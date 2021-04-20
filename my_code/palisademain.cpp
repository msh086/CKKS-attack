#include "palisade.h"

using namespace lbcrypto;

int main() {
    std::boolalpha(std::cout);

    // whether the sk holder uses insecure decoding algorithm. i.e image parts are kept and no extra noise is added
    bool insecure_client_encoding = true;

    // Step 1: Setup CryptoContext

    // A. Specify main parameters
    /* A1) Multiplicative depth:
     * The CKKS scheme we setup here will work for any computation
     * that has a multiplicative depth equal to 'multDepth'.
     * This is the maximum possible depth of a given multiplication,
     * but not the total number of multiplications supported by the
     * scheme.
     *
     * For example, computation f(x, y) = x^2 + x*y + y^2 + x + y has
     * a multiplicative depth of 1, but requires a total of 3 multiplications.
     * On the other hand, computation g(x_i) = x1*x2*x3*x4 can be implemented
     * either as a computation of multiplicative depth 3 as
     * g(x_i) = ((x1*x2)*x3)*x4, or as a computation of multiplicative depth 2
     * as g(x_i) = (x1*x2)*(x3*x4).
     *
     * For performance reasons, it's generally preferable to perform operations
     * in the shorted multiplicative depth possible.
     */
    uint32_t multDepth = 4;

    /* A2) Bit-length of scaling factor.
     * CKKS works for real numbers, but these numbers are encoded as integers.
     * For instance, real number m=0.01 is encoded as m'=round(m*D), where D is
     * a scheme parameter called scaling factor. Suppose D=1000, then m' is 10 (an
     * integer). Say the result of a computation based on m' is 130, then at
     * decryption, the scaling factor is removed so the user is presented with
     * the real number result of 0.13.
     *
     * Parameter 'scaleFactorBits' determines the bit-length of the scaling
     * factor D, but not the scaling factor itself. The latter is implementation
     * specific, and it may also vary between ciphertexts in certain versions of
     * CKKS (e.g., in EXACTRESCALE).
     *
     * Choosing 'scaleFactorBits' depends on the desired accuracy of the
     * computation, as well as the remaining parameters like multDepth or security
     * standard. This is because the remaining parameters determine how much noise
     * will be incurred during the computation (remember CKKS is an approximate
     * scheme that incurs small amounts of noise with every operation). The
     * scaling factor should be large enough to both accommodate this noise and
     * support results that match the desired accuracy.
     */
    uint32_t scaleFactorBits = 40;

    /* A3) Number of plaintext slots used in the ciphertext.
     * CKKS packs multiple plaintext values in each ciphertext.
     * The maximum number of slots depends on a security parameter called ring
     * dimension. In this instance, we don't specify the ring dimension directly,
     * but let the library choose it for us, based on the security level we
     * choose, the multiplicative depth we want to support, and the scaling factor
     * size.
     *
     * Please use method GetRingDimension() to find out the exact ring dimension
     * being used for these parameters. Give ring dimension N, the maximum batch
     * size is N/2, because of the way CKKS works.
     */
    uint32_t batchSize = 1 << 15;

    /* A4) Desired security level based on FHE standards.
     * This parameter can take four values. Three of the possible values
     * correspond to 128-bit, 192-bit, and 256-bit security, and the fourth value
     * corresponds to "NotSet", which means that the user is responsible for
     * choosing security parameters. Naturally, "NotSet" should be used only in
     * non-production environments, or by experts who understand the security
     * implications of their choices.
     *
     * If a given security level is selected, the library will consult the current
     * security parameter tables defined by the FHE standards consortium
     * (https://homomorphicencryption.org/introduction/) to automatically
     * select the security parameters. Please see "TABLES of RECOMMENDED
     * PARAMETERS" in  the following reference for more details:
     * http://homomorphicencryption.org/wp-content/uploads/2018/11/HomomorphicEncryptionStandardv1.1.pdf
     */
    SecurityLevel securityLevel = HEStd_NotSet;

    // The following call creates a CKKS crypto context based on the
    // arguments defined above.
    CryptoContext<DCRTPoly> cc =
            CryptoContextFactory<DCRTPoly>::genCryptoContextCKKS(
                    multDepth, scaleFactorBits, batchSize, securityLevel, 1 << 16);

    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension()
              << std::endl
              << std::endl;

    // Enable the features that you wish to use
    cc->Enable(ENCRYPTION);
    cc->Enable(SHE);

    // B. Step 2: Key Generation
    /* B1) Generate encryption keys.
     * These are used for encryption/decryption, as well as in generating
     * different kinds of keys.
     */
    /*
     * ckks.cpp line 33
     * */
    auto keys = cc->KeyGen();

    /* B2) Generate the relinearization key
     * In CKKS, whenever someone multiplies two ciphertexts encrypted with key s,
     * we get a result with some components that are valid under key s, and
     * with an additional component that's valid under key s^2.
     *
     * In most cases, we want to perform relinearization of the multiplicaiton
     * result, i.e., we want to transform the s^2 component of the ciphertext so
     * it becomes valid under original key s. To do so, we need to create what we
     * call a relinearization key with the following line.
     */
    cc->EvalMultKeyGen(keys.secretKey);

    /* B3) Generate the rotation keys
     * CKKS supports rotating the contents of a packed ciphertext, but to do so,
     * we need to create what we call a rotation key. This is done with the
     * following call, which takes as input a vector with indices that correspond
     * to the rotation offset we want to support. Negative indices correspond to
     * right shift and positive to left shift. Look at the output of this demo for
     * an illustration of this.
     *
     * Keep in mind that rotations work on the entire ring dimension, not the
     * specified batch size. This means that, if ring dimension is 8 and batch
     * size is 4, then an input (1,2,3,4,0,0,0,0) rotated by 2 will become
     * (3,4,0,0,0,0,1,2) and not (3,4,1,2,0,0,0,0). Also, as someone can observe
     * in the output of this demo, since CKKS is approximate, zeros are not exact
     * - they're just very small numbers.
     */
    cc->EvalAtIndexKeyGen(keys.secretKey, {1, -2});

    // Step 3: Encoding and encryption of inputs

    // Inputs
    vector<double> x1 = {0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0, 5.0};
    vector<double> x2 = {5.0, 4.0, 3.0, 2.0, 1.0, 0.75, 0.5, 0.25};

    // Encoding as plaintexts
    /*
     * cryptocontext.h line 1805 CryptoContextImpl::MakeCKKSPackedPlaintext
     * cryptocontext.h line 1785 CryptoContextImpl::MakeCKKSPackedPlaintext
     * ckkspackedencoding.cpp line 113 CKKSPackedEncoding::Encode
     * NOTE: the image part of complex values in the input vector are cleared out before encoding
     *  consequently we need to implement an insecure version of encoding
     *  due to polymorphism, we cannot add a new method or add a new parameter
     *  maybe use member fields
     * */
    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x1);
    Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(x2);

    std::cout << "Input x1: " << ptxt1 << std::endl;
    std::cout << "Input x2: " << ptxt2 << std::endl;

    // Encrypt the encoded vectors
    /*
     * ckks-impl.cpp line 716 LPAlgorithmCKKS<DCRTPoly>::Encrypt
     * */
    auto c1 = cc->Encrypt(keys.publicKey, ptxt1);
    auto c2 = cc->Encrypt(keys.publicKey, ptxt2);

    // Step 4: Evaluation

    // Homomorphic addition
    /*
     * NOTE: DCRTPoly := DCRTPolyImpl<BigVector>
     * ckks.cpp line 271 LPAlgorithmSHECKKS<Element>::EvalAdd // with Element = DCRTPoly
     * ckks.cpp line 147 LPAlgorithmSHECKKS<Element>::EvalAddCore // with Element = DCRTPoly
     * dcrtpoly.h line 1446 operator+ on DCRTPolyType := DCRTPolyImpl<VecType> // with VecType = BigVector
     * dcrtpoly.cpp line 626 DCRTPolyImpl<VecType>::Plus // with VecType = BigVector
     * NOTE: the data of DCRTPoly is held in an vector of PolyType := PolyImpl<NativeVector>
     * poly.h line 787 operator+ on PolyImpl<VecType> // with VecType = NativeVector
     * poly.cpp line 576 PolyImpl<VecType>::Plus // with VecType = NativeVector
     * mubintvecnat.cpp line 258 NativeVector<IntegerType>::ModAdd // with IntegerType = NativeInteger64
     * NOTE: NativeVector := bigintnat::NativeVector<NativeInteger64>
     *  NativeInteger64 := NativeIntegerT<uint64_t>
     *
     * NOTE: hierarchical view of the classes
     *  NativeInteger<T> provides (sometimes double-precision) arithmetic on native integer type T (here it is uint64_t)
     *  NativeVector<NativeInteger<T>> provides element-wise arithmetic on array of native integers w.r.t some modulus
     *  PolyImpl<NativeVector<NativeInteger<T>>> provides arithmetic on a single tower(NTT, sampling, etc.)
     *  DCRTPolyImpl<BigVector> provides arithmetic on double CRT ring elements(CRT & inv-CRT, multi-precibuiltsion, etc.)
     *
     * NOTE: built-in inverse computation
     *  DCRTPolyImpl<VecType>::MultiplicativeInverse
     * */
    auto cAdd = cc->EvalAdd(c1, c2);
    // Homomorphic subtraction
    auto cSub = cc->EvalSub(c1, c2);

    // Homomorphic scalar multiplication
    auto cScalar = cc->EvalMult(c1, 4.0);

    // Homomorphic multiplication, relinearization is done automatically
    auto cMul = cc->EvalMult(c1, c2);

    // Homomorphic rotations
    auto cRot1 = cc->EvalAtIndex(c1, 1);
    auto cRot2 = cc->EvalAtIndex(c1, -2);

    // Step 5: Decryption and output
    Plaintext result;
    // We set the cout precision to 8 decimal digits for a nicer output.
    // If you want to see the error/noise introduced by CKKS, bump it up
    // to 15 and it should become visible.
    std::cout.precision(8);
    std::cout << std::endl
              << "Results of homomorphic computations: " << std::endl;

    // Decrypt the result of addition
    /*
     * cryptocontext-impl.cpp line 43 CryptoContextImpl<DCRTPoly>::Decrypt
     *  decryption at line 71
     *      ckks-impl.cpp line 960 LPAlgorithmCKKS<DCRTPoly>::Decrypt
     *  decoding at line 90
     *      ckkspackedencoding.cpp line 330 CKKSPackedEncoding::Decode
     *
     * Decoding:
     *  defense is already present
     *  the stddev in Gaussian Distribution is 2 * (the stddev of the image part in message vector)
     *  the image part is cleared at the end of decoding for defense
     *
     * NOTE: defense noise
     *  defense noise is added in the following manner
     *  suppose decoded vector z corresponds to the ring element msg = (m_0, m_1, ..., m_(n-1))
     *  then conj(z) corresponds to msg' = (m_0, -m_(n-1), ..., -m_1)
     *  2 * Image(z) = z - conj(z) corresponds to msg - msg' = (0, m_1 + m_(n-1), m_2 + m_(n-2), ..., m_1 + m_(n-1))
     *  since only REAL vectors are accepted as inputs, the image part of z can somehow represent the error
     *  let the stddev be the standard deviation of the n-1 non-zero coefficients in (msg - msg')/2
     *  then the standard deviation for gaussian sampling is 2 * stddev
     *  besides, 2 * Real(z) = z + conj(z) corresponds to msg + msg', which is used as the actual input for decoding
     *  first (msg + msg')/2 is computed, then gaussian noise is added to each coefficient independently
     *
     * Full description of the defense
     * (the polynomial coefficients are already in multi-precision representation)
     * 1. the coefficients are rescaled s.t. the rescaling factor is exactly 2^p (where p is given by the user)
     * 2. compute the stddev of the n-1 non-zero coefficients of (msg - msg')/2
     * 3. if stddev < sqrt(n)/8, stddev := sqrt(n)/8 FIXME: why?
     * 4. if log2(stddev) > p - 10, throw an exception that precision is less than 10 bits FIXME: why?
     * 5. stddev := (M + 1) * stddev, where M defaults to 1
     * 6. compute (msg + msg')/2
     * 7. add gaussian noise to each coefficient in the polynomial above
     * 8. perform decoding fft (canonical embedding)
     * 9. clear the image part of the decoded vector
     * 10. estimate the approximate error as round(log2(stddev * sqrt(n)))
     * */
    cc->Decrypt(keys.secretKey, cAdd, &result);
    result->SetLength(batchSize);
    /*
     * some unrelated notes on the stream operator
     * 1. the standard library provides us with a templated stream operator on shared_ptr
     *  its signature is
     *      template <class T, class U, class V>
     *      std::basic_ostream<U, V>& operator<<(std::basic_ostream<U, V>& os, const std::shared_ptr<T>& ptr);
     * 2. the user-defined streams operator has a signature of
     *      std::ostream& operator<<(std::ostream& out, const Plaintext item)
     * 3. both of them are in scope (and the former can be found by ADL)
     *  but the second one is more specified in that ostream is a specification of basic_ostream
     *  so the user-defined version is selected by the compiler
     * */
    std::cout << "Estimated precision in bits: " << result->GetLogPrecision()
              << std::endl;

    // Decrypt the result of subtraction
    cc->Decrypt(keys.secretKey, cSub, &result);
    result->SetLength(batchSize);

    // Decrypt the result of scalar multiplication
    cc->Decrypt(keys.secretKey, cScalar, &result);
    result->SetLength(batchSize);

    // Decrypt the result of multiplication
    cc->Decrypt(keys.secretKey, cMul, &result);
    result->SetLength(batchSize);

    // Decrypt the result of rotations
    cc->Decrypt(keys.secretKey, cRot1, &result);
    result->SetLength(batchSize);
    std::cout
            << std::endl
            << "In rotations, very small outputs (~10^-10 here) correspond to 0's:"
            << std::endl;

    auto &victim = cRot2;

    if (insecure_client_encoding)
        cc->set_insecure_encodings(1);
    cc->Decrypt(keys.secretKey, victim, &result);
    result->SetLength(batchSize);


    // attack
    // decrypted polynomial
    Plaintext decrypted = cc->GetPlaintextForDecrypt(
            victim->GetEncodingType(), victim->GetElements()[0].GetParams(),
            cc->GetEncodingParams());
    cc->GetEncryptionAlgorithm()->Decrypt(keys.secretKey, victim, &decrypted->GetElement<Poly>());
    const auto &msg = decrypted->GetElement<Poly>();

    // re-encoding
    cc->set_insecure_encodings(1);
    auto re_encoded = cc->MakeCKKSPackedPlaintext(result->GetCKKSPackedValue(), result->GetDepth(), result->GetLevel(),
                                                  victim->GetCryptoParameters()->GetElementParams());
    auto &sk_ele = keys.secretKey->GetPrivateElement();
    auto &re_encoded_ele = re_encoded->GetElement<DCRTPoly>();
    auto re_encoded_ele_big = re_encoded_ele.CRTInterpolate();

    std::cout << "re-encoding success? " << (msg == re_encoded_ele_big) << std::endl;
    std::cout << "max encoding error= " << (msg - re_encoded_ele_big).Norm() << std::endl;

    // key recovery
    auto &m = re_encoded_ele;
    auto b = victim->GetElements()[0], a = victim->GetElements()[1];
    m -= b;
    m *= a.MultiplicativeInverse();
    auto sk_copy = sk_ele;
    sk_copy.DropLastElements(sk_copy.GetNumOfElements() - m.GetNumOfElements());
    std::cout << "recovery ok? " << (sk_copy == m) << std::endl;
    return 0;
}
