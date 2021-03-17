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
#include <ctime>

void printParams(const Plaintext& plaintext) {
    printf("log(p) = %ld, log(q) = %ld\n", plaintext.logp, plaintext.logq);
}

/**
 * max([valid bits of coeff for coeff in arr[0:n]])
 * */
long maxBits(const NTL::ZZ* arr, long n) {
    long res = 0;
    for(long i = 0; i < n; i++){
        res = std::max(res, NTL::NumBits(arr[i]));
    }
    return res;
}

void printPlainMag(const Plaintext& plaintext) {
    printf("magnitude of %ld\n", maxBits(plaintext.mx, plaintext.n));
}

/**
 * convert value in [0, 2^logQ - 1] to [-2^(logQ - 1), 2^(logQ - 1) - 1]
 * */
void twosComplementToSigned(Plaintext& plaintext) {
    NTL::ZZ modulus = NTL::ZZ(1) << plaintext.logq;
    for(long i = 0; i < plaintext.n; i++) {
        if (NTL::NumBits(plaintext.mx[i]) == plaintext.logq)
            plaintext.mx[i] -= modulus;
    }
}

std::complex<double> sampleUnitCircle() {
    unsigned long bits = NTL::RandomWord(); // this only works on 64-bit machines
    double squareR = (double)bits / UINT64_MAX;
    bits = NTL::RandomWord();
    double theta = (double)bits / UINT64_MAX * 2 * M_PI;
    return std::polar(std::sqrt(squareR), theta);
}

std::complex<double>* sampleUnitCircleArr(int n) {
     auto *vec = new std::complex<double>[n]; // bad habit, but conforms with HEAAN :)
     for(int i = 0; i < n; i++)
         vec[i] = sampleUnitCircle();
     return vec;
}

NTL::ZZX composeZZArray(NTL::ZZ* arr, long n) {
    NTL::ZZX res;
    res.SetLength(n + 1);
    for(long i = 0; i < n; i++)
        res[i] = arr[i];
    return res; // RVO will help us
}

NTL::ZZ_pE convRing(const NTL::ZZX& src) {
    return NTL::conv<ZZ_pE>(NTL::conv<ZZ_pX>(src));
}

bool checkSame(const Plaintext& p1, const Plaintext& p2) {
    long n = p1.n;
    if (n != p2.n)
        return false;
    for(long i = 0; i < n; i++) {
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
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
    scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message
    
    // Make Random Array of Complex //
    std::complex<double>* mvec1 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);
    std::complex<double>* mvec2 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);

    // logp seems to be the scaling factor, while logq is the ciphertext modulus
    // Encrypt Two Arry of Complex //
    Plaintext plain1 = scheme.encode(mvec1, slots, logp, logQ);
    printParams(plain1);
    printPlainMag(plain1);

    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);

    // Addition //
    Ciphertext cipherAdd = scheme.add(cipher1, cipher2);

    // Multiplication And Rescale //
    Ciphertext cipherMult = scheme.mult(cipher1, cipher2);
    Ciphertext cipherMultAfterReScale = scheme.reScaleBy(cipherMult, logp);

    // Rotation //
    long idx = 1;
    Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);

    // choose attack victim //
    Ciphertext& victim = cipher1;

    // Decrypt //
    Plaintext dmsg1 = scheme.decryptMsg(secretKey, victim);
    twosComplementToSigned(dmsg1);
    std::complex<double>* dvec1 = scheme.decode(dmsg1);

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

    auto m_ring = convRing(composeZZArray(re_encoded.mx, ring.N));
    auto b_ring = convRing(composeZZArray(victim.bx, ring.N));
    auto a_ring = convRing(composeZZArray(victim.ax, ring.N));

    m_ring -= b_ring;
    m_ring /= a_ring;

    auto s_ring = convRing(composeZZArray(secretKey.sx, ring.N));
    std::cout << "recovery ok ? " << bool(s_ring == m_ring) << std::endl;

    return 0;

}