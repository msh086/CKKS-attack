//
// Created by msh on 2021/3/10.
//

// HEAAN includes code like "using namespace std", which is really a bad habit
#include "HEAAN.h"
#include <iostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/GF2X.h> // ZZ_pX when p = 2
#include <NTL/GF2E.h> // ZZ_pE when p = 2
#include <ctime>
#include <exception>

#define PRINTERROR(err) std::cerr << "In function " << __FUNCTION__  << " at line " << __LINE__ << ": "\
    << e.what() << std::endl;


void printParams(const Plaintext &plaintext) {
    printf("log(p) = %ld, log(q) = %ld\n", plaintext.logp, plaintext.logq);
}

/**
 * max([valid bits of coeff for coeff in arr[0:n]])
 * */
long maxBits(const NTL::ZZ *arr, long n) {
    long res = 0;
    for (long i = 0; i < n; i++) {
        res = std::max(res, NTL::NumBits(arr[i]));
    }
    return res;
}

void printPlainMag(const Plaintext &plaintext) {
    printf("magnitude of %ld\n", maxBits(plaintext.mx, plaintext.n));
}

/**
 * convert value in [0, 2^logQ - 1] to [-2^(logQ - 1), 2^(logQ - 1) - 1]
 * */
void unsignedToSigned(Plaintext &plaintext) {
    NTL::ZZ modulus = NTL::ZZ(1) << plaintext.logq;
    for (long i = 0; i < plaintext.n; i++) {
        if (NTL::NumBits(plaintext.mx[i]) == plaintext.logq)
            plaintext.mx[i] -= modulus;
    }
}

/**
 * return a complex value sampled uniformly from the unit circle
 * */
std::complex<double> sampleUnitCircle() {
    unsigned long bits = NTL::RandomWord(); // this only works on 64-bit machines
    double squareR = (double) bits / UINT64_MAX;
    bits = NTL::RandomWord();
    double theta = (double) bits / UINT64_MAX * 2 * M_PI;
    return std::polar(std::sqrt(squareR), theta);
}

/**
 * return a complex-valued array, all elements sampled uniformly from the unit circle
 * */
std::complex<double> *sampleUnitCircleArr(int n) {
    auto *vec = new std::complex<double>[n]; // bad habit, but conforms with HEAAN :)
    for (int i = 0; i < n; i++)
        vec[i] = sampleUnitCircle();
    return vec;
}

/**
 * convert ZZ-array-styled polynomials used in HEAAN into ZZX
 * */
NTL::ZZX composeZZArray(NTL::ZZ *arr, long n) {
    NTL::ZZX res;
    res.SetLength(n + 1);
    for (long i = 0; i < n; i++)
        res[i] = arr[i];
    return res; // RVO will help us
}

/**
 * convert element from ZZX(Z[x]) into ZZ_pE(R_q)
 * */
NTL::ZZ_pE convRing(const NTL::ZZX &src) {
    return NTL::conv<ZZ_pE>(NTL::conv<ZZ_pX>(src));
}

/**
 * convert element from ZZ_pE(R_q) into ZZX(Z[x])
 * */
NTL::ZZX deconvRing(const NTL::ZZ_pE& src){
    return NTL::conv<ZZX>(NTL::conv<ZZ_pX>(src));
}

/**
 * convert element in R_q into R_2 (q is a power of 2)
 * */
NTL::GF2E convModDown(const NTL::ZZ_pE &src) {
    return NTL::conv<NTL::GF2E>(
            NTL::conv<NTL::GF2X>(
                    NTL::conv<NTL::ZZX>(
                            NTL::conv<NTL::ZZ_pX>(src)
                    )
            )
    );
}

/**
 * convert element in R_2 into R_q (q is a power of 2)
 * */
NTL::ZZ_pE convModUp(const NTL::GF2E &src) {
    return NTL::conv<NTL::ZZ_pE>(
            NTL::conv<NTL::ZZ_pX>(
                    NTL::conv<NTL::ZZX>(
                            NTL::conv<NTL::GF2X>(src)
                    )
            )
    );
}

/**
 * b + a * s = m, with a, b, s, m in R_q
 * recover s by computing a^-1 * (m - b)
 * NOTE: this requires the a is invertible in R_q
 * */
NTL::ZZ_pE recoverByInv(const NTL::ZZ_pE& m, const NTL::ZZ_pE& a, const NTL::ZZ_pE& b) {
    NTL::ZZ_pE res = m;
    res -= b;
    try {
        res /= a;
    } catch (NTL::ErrorObject& e) {
        PRINTERROR(e)
    }
    return res;
}

/**
 * b + a * s = m, with a, b, s, m in R_q
 * recover s using the same trick that Li and Micciancio have adopted
 * this trick takes advantage of the fact that sk_i = -1, 0, 1
 * NOTE: this requires the (a mod 2) is invertible in R_2
 * */
NTL::ZZ_pE recoverByTrick(const NTL::ZZ_pE& m, const NTL::ZZ_pE& a, const NTL::ZZ_pE& b) {
    auto m2 = convModDown(m), a2 = convModDown(a), b2 = convModDown(b);
    // compute [a]^-1
    NTL::GF2E a2inv;
    try {
        NTL::inv(a2inv, a2);
    } catch (NTL::ErrorObject& e) {
        PRINTERROR(e);
        return NTL::ZZ_pE();
    }
    // compute [s]
    auto s2 = (m2 - b2) * a2inv;
    // compute h, where 2 * h + [s] = s
    auto ah = convRing(deconvRing(m - b - a * convModUp(s2)) / 2);
    auto h2 = convModDown(ah) * a2inv;
    return (-2) * convModUp(h2) + convModUp(s2);
}


bool checkSame(const Plaintext &p1, const Plaintext &p2) {
    long n = p1.n;
    if (n != p2.n)
        return false;
    for (long i = 0; i < n; i++) {
        if (p1.mx[i] != p2.mx[i]) {
            std::cout << "difference at " << p1.mx[i] << " " << p2.mx[i] << std::endl;
            return false;
        }
    }
    return true;
}

int main() {
    NTL::ZZ seed(time(nullptr));
    std::cout << "seed is " << seed << std::endl;
    NTL::SetSeed(seed);

    // Parameters //
    long logN = 15;
    long logQ = 353;
    long logp = 30; ///< Larger logp will give you more correct result (smaller computation noise)
    long slots = 1 << (logN - 1); ///< This should be power of two
    long numThread = 8;

    // Construct and Generate Public Keys //
    TimeUtils timeutils;
    Ring ring(logN, logQ);
    /**
     * hamming weight defaults to 64
     * each non-zero coefficient of secret key has value of 1 or -1 with uniform possibility
     * */
    SecretKey secretKey(ring);

    int total = 1, inv_success = 0, trick_success = 0;

    for(int iter = 0; iter < total; iter++) {
        Scheme scheme(secretKey, ring);
        scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
        scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message

        // Make Random Array of Complex //
        std::complex<double> *mvec1 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);
//        std::complex<double> *mvec2 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);

        // logp seems to be the scaling factor, while logq is the ciphertext modulus
        // Encrypt Two Arry of Complex //
        Plaintext plain1 = scheme.encode(mvec1, slots, logp, logQ);
        printParams(plain1);
        printPlainMag(plain1);

        Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
//        Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);

        // Addition //
//        Ciphertext cipherAdd = scheme.add(cipher1, cipher2);

        // Multiplication And Rescale //
//        Ciphertext cipherMult = scheme.mult(cipher1, cipher2);
//        Ciphertext cipherMultAfterReScale = scheme.reScaleBy(cipherMult, logp);

        // Rotation //
//        long idx = 1;
//        Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);

        // choose attack victim //
        Ciphertext &victim = cipher1;

        // Decrypt //
        Plaintext dmsg1 = scheme.decryptMsg(secretKey, victim);
        unsignedToSigned(dmsg1);
        std::complex<double> *dvec1 = scheme.decode(dmsg1);

        // attack
        // HEAAN will scale the inv-embedded polynomial by 2^(logp + ring.logQ)
        // we want to get rid of the extra scaling factor 2^(ring.logQ)
        // although round(round(m * 2^(logp + ring.logQ)) / 2^ring.logQ) will be the same as round(m * 2^logp)
        // unless the fractional part of m * 2^logp has a leading 0-bit, followed by logQ 1-bits
        // i.e m * 2^logp mod 1 = b'.01...1xxxx', where the count of bit 1 is logQ (xxx means any)
        // here is a simple trick to achieve this, without modifying library code
        Plaintext re_encoded = scheme.encode(dvec1, slots, dmsg1.logp - ring.logQ, dmsg1.logq);
        re_encoded.logp += ring.logQ;

        printParams(dmsg1);
        printPlainMag(dmsg1);
        printPlainMag(re_encoded);

        std::cout << "re-encoding ok ? " << std::boolalpha << checkSame(re_encoded, dmsg1) << std::endl;

        // init rings
        NTL::ZZ_p::init(ring.Q);
        NTL::ZZX PhiM;
        PhiM.SetLength(ring.N + 1);
        PhiM[0] = 1;
        PhiM[ring.N] = 1;
        NTL::ZZ_pE::init(NTL::conv<NTL::ZZ_pX>(PhiM));
        NTL::GF2E::init(NTL::conv<NTL::GF2X>(PhiM));

        auto m_ring = convRing(composeZZArray(re_encoded.mx, ring.N));
        auto b_ring = convRing(composeZZArray(victim.bx, ring.N));
        auto a_ring = convRing(composeZZArray(victim.ax, ring.N));

        auto recover = recoverByInv(m_ring, a_ring, b_ring);

        auto s_ring = convRing(composeZZArray(secretKey.sx, ring.N));
        std::cout << "recovery by inverse ok ? " << bool(s_ring == recover) << std::endl;

        inv_success += s_ring == recover;
        /**
         * maybe a is non-invertible, but when all of sk's coefficients are constrained into {-1, 0, 1}
         * we can found a way to decide such sk without computing a^-1
         * the method below is taken from CKKSKeyRecovery by Li and Micciancio
         *
         * 1. we know b + a * s = m, since ciphertext modulus q = 2^k,
         *  natural reduction f: R_q => R_2 forms a homomorphism (although a single-way one)
         *  ### Note: we will note the image of a polynomial 'm' in R_q under such mapping '[m]'
         *  so the following equation holds:
         *  [b] + [a] * [s] = [m]   (Eq.1)
         *  FIXME the statement below can be more precise, see the discussion at the end of this comment block
         *  the ring R_2 is friendlier than R_q because the set of invertible elements in R_q is a subset of
         *  the set of invertible elements in R_2
         *  so with great possibility we can compute [s] = ([m] - [b]) * [a]^-1
         *  ### Note: for a polynomial 'm', we index its coefficient w.r.t x^i with 'm_i'
         *  since sk's coefficients compose of {-1, 0, 1}, [s]_i = 0 iff sk[i] = 0, and [s]_i = 1 iff sk[i] = 1 or -1
         *
         * 2. we have already distinguished all non-zero coefficients in sk, now it only remains to
         *  decide whether a non-zero coefficient is -1 or 1
         *  ### Note: we will not distinguish between [m] in R_2 and [m] in R_q, which meaning it takes depends on
         *  ###     with whom [m] is multiplied / added. e.g. in [m] * n, [m] is in R_q, while in [m] * [n], [m] is in R_2
         *  subtracting [s] from s, we have
         *  b + a * (h + [s]) = m   (Eq.2)
         *  where only h is unknown and is what we want
         *  (alternatively, we can add [s] to s, like Li and Micciancio did, which is similar)
         *  we know h_i = 0 iff s_i = 1 or 0, h_i = -2 (or 2^k - 2) iff s_i = -1
         *  all coefficients of h are even, so we can divide it by 2
         *  hx := h / 2
         *  which means
         *  b + a * (2 * hx + [s]) = m  ==>>
         *  a * hx = (m - b - a * [s]) / 2  (Eq.3)
         *  now hx_i = 0 iff s_i = 1 or 0, hx_i = -1 (or 2^(k - 1) - 1) iff s_i = -1
         *  mapping Eq.3 into R_2, we get
         *  [a] * [hx] = [(m - b - a * [s]) / 2]
         *  from which we can easily compute [hx] by multiplying both sides with [a]^-1
         *
         * 3. sk = -2 * [hx] + [s]
         *  if we defined h := s + [s], then here sk = 2 * [hx] - [s]
         *
         * 4.
         * FIXME strict proof needed
         *   a invertible polynomial 'm' in R_2 must satisfy that m(1) = 1
         *   since for any 'm' in R_2, f: R_2 -> GF(2), f(m) = m(1) is a homomorphic mapping
         *   (because 1 is a root of modulus polynomial x^n + 1)
         *   suppose m is invertible, i.e. m * m^-1 = 1, we apply f to both sides,
         *   we get m(1) * m^-1(1) = 1, then m(1) = m^-1(1) = 1
         *   what remains to be proved is that set S = {m in R_2, m(1) = 1 and m is non-invertible} is empty
         *   below gives a not-so-strict proof on this
         * Proof(sketch):
         *  consider the set S mentioned above
         *  a polynomial 'm' in R_2 is non-invertible iff gcd(x^n + 1, m) != 1
         *  while x^n + 1 = (x + 1)^n in GF(2), this can be proved by induction
         *  so if m is non-invertible, x + 1 | m
         *  we write m = (x + 1) * a, where a is some polynomial in R_2
         *  from the discussion above, we know
         *  m(1) = (x + 1)(1) * a(1) =  (1 + 1)(1) * a(1) = 0
         *  this contradicts with the assumption that m(1) = 1
         *  so the set S is empty
         *  which indicates that m is invertible iff m(1) = 1
         *  FIXME: I'm not sure if every step is correct in GF(2), for example
         *   is factorizing a polynomial in GF(2)[x] the same as that in Z[x]?
         * */
        auto recoverTrick = recoverByTrick(m_ring, a_ring, b_ring);
        std::cout << "recovery by trick ok ? " << bool(s_ring == recoverTrick) << std::endl;

        trick_success += s_ring == recoverTrick;

        delete[] mvec1;
    }

    printf("inv success: %d\ntrick success: %d\ntotal: %d\n", inv_success, trick_success, total);

    return 0;

}