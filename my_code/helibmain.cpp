//
// Created by msh on 2021/1/27.
//
#include <iostream>
#include <vector>
#include <helib/helib.h>
#include <NTL/ZZ_pE.h>

using namespace std;
using namespace helib;

/**
 * decode a polynomial into a complex vector
 * this function is copied from EaCx.cpp::CKKS_decode line 33
 * it is declared as a global static function so we have to reimplement it here
 * */
void CKKS_decode(const NTL::ZZX& pp,
                 NTL::xdouble xfactor,
                 const PAlgebra& palg,
                 std::vector<cx_double>& ptxt)
{
    const long MAX_BITS = 400;
    long nBits = NTL::MaxBits(pp) - MAX_BITS;
    double factor;

    // This logic prevents floating point overflow
    if (nBits <= 0) {
        CKKS_canonicalEmbedding(ptxt, pp, palg);
        factor = NTL::to_double(xfactor);
    } else {
        long dpp = deg(pp);
        std::vector<double> pp_scaled(dpp + 1);
        NTL::ZZ tmp;
        for (long i : range(dpp + 1)) {
            RightShift(tmp, pp.rep[i], nBits);
            pp_scaled[i] = NTL::to_double(tmp);
        }
        CKKS_canonicalEmbedding(ptxt, pp_scaled, palg);
        factor = NTL::to_double(xfactor / NTL::power2_xdouble(nBits));
    }

    for (cx_double& cx : ptxt) // divide by the factor
        cx /= factor;
}

int main(int argc, char** argv){

    Context context = ContextBuilder<CKKS>()
            .m(1 << 17) // n = 2^12
            .bits(1098) // settings for bits and c0 follow recommendation in CKKS tutorial
            .c(2)
            .precision(20)
            .build();
    cout << "security level = " << context.securityLevel() << "\n";

    long n = context.getNSlots();

    SecKey secretKey(context);

    secretKey.GenSecKey();

    // NOTE: reference initialization: publicKey is essentially secretKey
    //  encryption will actually use secret key
    //  to avoid such situation, omit & to use copy construction
    const PubKey& publibKey = secretKey;

    vector<double> vec(n);
    for(size_t i = 0;i < n ;i++)
        vec[i] = sin(2 * PI * i / n);

    /**
     * plaintext is not encoded in PtxtArray
     * it simply holds a the unencoded vector
     * */
    PtxtArray ptxt0(context, vec);
    Ctxt c0(publibKey);
    /**
     * Encoding is performed here in PtxtArray::encrypt
     * the underlying encoding function is norms.cpp::CKKS_embedInSlots line 574-615
     * the coefficients of the encoded polynomial are stored in a NTL::vec<long> object
     *
     * calling sequence:
     * EncryptedArray.h::EncryptedArray::encode line 1727 -> call on inner generic pointer
     * EaCx.cpp::EncryptedArrayDerived<PA_cx>::encode line 280 -> get complex vector in PlaintextArray
     * EaCx.cpp::EncryptedArrayDerived<PA_cx>::encode line 238 -> determine scaling factor
     * norms.cpp::CKKS_embedInSlots line 574 -> actual canonical embedding
     *
     * the scaling factor is automatically chosen
     *  1. assume the encoding error (for each polynomial coefficient) is uniform in [-0.5, 0.5]
     *  2. assume each slot after decoding follows a normal distribution with the same variance
     *     here the variance is N*sigma^2, where sigma^2 is the variance for each polynomial coefficient
     *     which would be Q^2/12 in the case of uniform continuous / discrete distribution in range [0, q]
     *  3. the probability that a uniform random variable with variance sigma^2
     *     falls out of the range of [-x, x] is erfc(x/(sigma*sqrt(2))).
     *     then the probability that it falls out of [-scale*x, scale*x] is delta = erfc(x/sqrt(2))
     *     (here the scale has nothing to do with the scaling factor).
     *     given a vector of phi(m) elements, its infinity norm is within scale*x
     *     is no more than delta * phi(m).
     *     so we can get a high probability bound on the L-infty norm of decoded noise (slot error)
     *  4. choose a scaling factor s.t. the error on each slot is not more than 2^(-precision)
     *   i.e. slot_err / scale <= 2^(-precision)
     * NOTE: in current implementation, slot_err = 10 * sqrt(N/12)
     *  let N = 2^k
     *  scale = 2^ceil(log_2(2^precision * slot_err)) = 2^ceil(precision + k/2 + log_2(5/sqrt(3)))
     *  where log_2(5/sqrt(3)) is approximately 1.5
     * */
    ptxt0.encrypt(c0);

    PtxtArray ptxt1(context);
    ptxt1.random();
    Ctxt c1(publibKey);
    ptxt1.encrypt(c1);

    Ctxt c2 = c0;
    c2 += c1;

    PtxtArray out(context);
    /**
     * Extra noise is added to plaintext polynomial after decryption to mitigate CKKS attack
     * see EaCx.cpp::EncryptedArray<PA_cx>::decrypt line 88 and Ctxt.cpp::addedNoiseForCKKSDecryption line 3026
     * the PRNG seed is hashed from the secret key and ciphertext
     *
     * in Ctxt.cpp::addedNoiseForCKKSDecryption, sigma is the standard deviation of each polynomial coefficient
     * sampleGaussianBoundedEffectiveBound returns a bound factor B s.t. the decoded vector is bounded by B*sigma
     * this function returns sqrt(N*ln(N))
     * TODO: to be honest I don't understand its logic, but I guess sqrt(ln(N)) plays the role of scale in erfc?
     *
     * here we need to decrypt without adding additional noise to perform attack
     * following code is similar to those in EaCx.cpp::EncryptedArray<PA_cx>::decrypt line 88
     * */
//    out.decrypt(c2, secretKey);
    vector<cx_double> vec_out;
//    out.store(vec_out);

    NTL::ZZX pp;
    secretKey.Decrypt(pp, c2);
    NTL::xdouble factor = c2.getRatFactor();
    auto palgebra = context.getAlMod().getZMStar();
    CKKS_decode(pp, factor, palgebra, vec_out);

    cout << "max bit length: " << NTL::NumBits(largestCoeff(pp)) << endl;
    cout << "scaling bits: " << NTL::NumBits(NTL::conv<NTL::ZZ>(factor)) << endl;

    out.load(vec_out);
    // re-encode
    zzX re_encoded;
    NTL::ZZX re_encoded_ZZX;
    CKKS_embedInSlots(re_encoded_ZZX, vec_out, palgebra, NTL::conv<double>(factor));
//    convert(re_encoded_ZZX, re_encoded);

    cout << "re-encoding ok? " <<(re_encoded_ZZX == pp ? "true" : "false") << endl;

    /// start calculation
    auto modulus = c2.getContext().productOfPrimes(c2.getContext().getCtxtPrimes());
    auto phimx = context.getZMStar().getPhimX();
    NTL::ZZ_p::init(modulus);
    NTL::ZZ_pX phim_mod_q;
    NTL::conv(phim_mod_q, phimx);

    NTL::ZZ_pE::init(phim_mod_q);
    NTL::ZZ_pE rhs;
    NTL::conv(rhs, NTL::conv<NTL::ZZ_pX>(re_encoded_ZZX));

    NTL::ZZX tmp;
    c2[0].toPoly(tmp);

    rhs -= NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(tmp));

    c2[1].toPoly(tmp);
    rhs /= NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(tmp));

    secretKey.getSK().toPoly(tmp);
    auto real_sk = NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(tmp));

    cout << "recovered ok? " << (rhs == real_sk ? "true" : "false") << endl;

//    PtxtArray plain_res = ptxt0;
//    plain_res *= ptxt1;
//
//    double distance = Distance(plain_res, out);
//    printf("distance = %e\n", distance);

    return 0;
}