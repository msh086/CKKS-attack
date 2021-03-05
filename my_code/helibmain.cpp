//
// Created by msh on 2021/1/27.
// NOTE: 2021/3/4
//  the attack on HElib doesn't work for this version (v2.0.0)
//  the scaling factor will grow so large, so that casting plaintext polynomial coefficients
//  into double data type will result in rounding error > 0.5
//  such type conversion is necessary for encoding / decoding
//
#include <iostream>
#include <vector>
#include <helib/helib.h>
#include <NTL/ZZ_pE.h>
#include "helibattack.h"
#include "utils.h"

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
    helibattack::print_NTL_macros();

    Context context = ContextBuilder<CKKS>()
            .m(1 << 17) // n = 2^16
            .bits(300) // settings for bits and c0 follow recommendation in CKKS tutorial
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

    vector<cx_double> vec(n);
    utils::sample_complex_circle(vec, n, 1);

    /**
     * plaintext is not encoded in PtxtArray
     * it simply holds the unencoded vector
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
     *     then the probability that it falls out of [-scale*sigma, scale*sigma] is delta = erfc(scale/sqrt(2))
     *     (here the scale has nothing to do with the scaling factor).
     *     given a vector of phi(m) elements, its infinity norm is within scale*x
     *     is no more than delta * phi(m).
     *     so we can get a high probability bound on the L-infty norm of decoded noise (slot error)
     *  4. choose a scaling factor s.t. the error on each slot is not more than 2^(-precision)
     *   i.e. slot_err / scale <= 2^(-precision)
     * NOTE: in current implementation, slot_err = 10 * sqrt(N/12), i.e. scale = 10, -log2(erfc(10/sqrt(2)) = 75.8
     *  let N = 2^k
     *  scale = 2^ceil(log_2(2^precision * slot_err)) = 2^ceil(precision + k/2 + log_2(5/sqrt(3)))
     *  where log_2(5/sqrt(3)) is approximately 1.5
     *
     * Encryption calling sequence depends on whether the PublicKey instance is actually a SecKey instance
     * in the case of true, keys.cpp line 1620 SecKey::Encrypt is called:
     *  a NEW RLWE instance is generated, playing the role of r <- ZO(0.5), e' <- GD(3.2), r * pk + e'
     * otherwise, keys.cpp line 756 PubKey::Encrypt is called:
     *  r <- ZO(0.5), e' <- GD(3.2), r * pk + e'
     * then in both functions, the error in the RLWE instance is estimated as 'err_bound'
     * while the encoding error is called 'err'
     * like that in homomorphic multiplication, HElib wants the precision loss in encryption is not more than 1 bit
     * i.e err_bound is no greater than err
     * so an extra scaling factor 'ef' is defined as ceil(err_bound / err)
     * then the plaintext polynomial is scaled up by ef, s.t. encryption noise is smaller than encoding noise
     * */
    ptxt0.encrypt(c0);

    PtxtArray ptxt1(context);
    ptxt1.random();
    Ctxt c1(publibKey);
    ptxt1.encrypt(c1);

    c1 *= c0;
    c1.dropSmallAndSpecialPrimes();

    helibattack::print_xdouble(c1.getRatFactor());
    /**
     * >> An attempt to analyze weird scaling up/down behaviour <<
     * in Ctxt.cpp::modSwitchAddedNoiseBound line 2535, this function calculates B_scale:
     *  c = (b, a) -->rescaling down by q--> c' = (round(b/q), round(a/q))
     *  e_scale = <c', sk> - <c, sk> / q = (b/q - round(b/q)) + (a/q - round(a/q)) * s
     *  here, we assume the rounding errors (b/q - round(b/q)) and (a/q - round(a/q))
     *  both follow a uniform distribution on [-0.5, 0.5]
     *  we note those rounding errors as random variable r_1, r_2 ~ U(-0.5, 0.5)
     *  according to inequalities on canonical norm, the infinite norm of plaintext error is
     *  B_scale = B(r_1 + r_2 * s) <= B(r_1) + B(r_2) * B(s) = (1 + B(s)) * noiseBoundUniform(0.5)
     *  // we are familiar with noiseBoundUniform, it returns a high-probability bound
     *  // on plaintext error, according to erfc and central limit theorem
     *
     * in Ctxt.cpp::computeIntervalForMul line 1589, look at document 5.3.2 for detailed explanation
     *  after rescaling from modulus q_i to q, the error bound will be q/q_i * B + B_scale
     *  to prevent error bound from growing too fast, HElib wants rescaled error (q/q_i * B) to
     *  be no greater than rounding error (B_scale), in this way the precision loss is only 1 bit per rescaling
     *  HElib decides q according to the rule above, so we have:
     *      log(q) \approx target = log(B_scale) + log(q_i) - log(B)
     *  what this function do is exactly calculating such log(q)
     *  Extra note: the document says HElib searches log(q) in [target - 4, target - 1]
     *  this is the case for BGV but not for CKKS; in the case of CKKS, the range is [target + 1, target + 4]
     *
     * in Ctxt.cpp::reLinearize line 706:
     *  first all small primes and special primes are dropped, this may affect ratFactor
     *  the special primes are added to prime set, causing the ratFactor to expand by the same ratio
     *  the special primes are NOT dropped after relinearization
     * */

    Ctxt c2 = c0;
    helibattack::eval_polynomial_inplace(utils::func_map.at(utils::F_EXPONENT), c2);
    utils::element_wise_eval_polynomial_inplace(utils::func_map.at(utils::F_EXPONENT), vec);

    PtxtArray out(context);
    /**
     * Extra noise is added to plaintext polynomial after decryption to mitigate CKKS attack
     * see EaCx.cpp::EncryptedArray<PA_cx>::decrypt line 88 and Ctxt.cpp::addedNoiseForCKKSDecryption line 3026
     * the PRNG seed is hashed from the secret key and ciphertext
     * the default std-dev is 6.4 (they multiplied 3.2 by 2 to get better security)
     *
     * in Ctxt.cpp::addedNoiseForCKKSDecryption, sigma is the standard deviation of each polynomial coefficient
     * sampleGaussianBoundedEffectiveBound returns a bound factor B s.t. the decoded vector is bounded by B*sigma
     * this function returns sqrt(N*ln(N))
     * FIXME complex gaussian tail bound probability needs proof
     * NOTE: see documentation 2.6.1
     *  each element of the decoded vector satisfies a complex gaussian distribution of variance sigma^2
     *  let such random variable be x
     *  then P(|x| > B) = exp(-B^2/sigma^2) /// need proof!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     *  let B = sigma * sqrt(log(n)), then P(|x| > B) = 1/n
     *  the heading half of the decoded vector is conjugate with the tailing half
     *  so we only need to consider n/2 elements
     *  the probability that the whole n/2-lengthed vector exceeds bound B is no greater than
     *  P(|x| > B) * n/2 = 1/2
     *  which means it's likely that the loop:
     *      sample a <- R_q
     *      decode v <- Dcd(a)
     *      test infty_norm(v) > B ?
     *          true -> go back
     *          false -> return
     *  can find a wanted a in no too much loops
     *  P.S. the limit on max loops is 1000, which means a runtime error caused by bounded gaussian sample failure
     *  is no more than 2^-1000
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

    cout << "max poly coeff bit length: " << NTL::NumBits(largestCoeff(pp)) << endl;
    cout << "scaling bits: " << NTL::NumBits(NTL::conv<NTL::ZZ>(factor)) << endl;

    cout << "max error: " << utils::max_error(vec, vec_out) << endl;
    cout << "rel error: " << utils::relative_error(vec, vec_out) << endl;

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