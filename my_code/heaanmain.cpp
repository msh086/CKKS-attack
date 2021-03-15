//
// Created by msh on 2021/3/10.
//

#include "HEAAN.h"

int main() {

    // Parameters //
    long logN = 15;
    long logQ = 353;
    long logp = 30; ///< Larger logp will give you more correct result (smaller computation noise)
    long slots = 1024; ///< This should be power of two
    long numThread = 8;

    // Construct and Generate Public Keys //
    TimeUtils timeutils;
    Ring ring(logN, logQ);
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
    scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message

    // Make Random Array of Complex //
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);

    // Encrypt Two Arry of Complex //
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

    // Decrypt //
    complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);
    complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);

    return 0;

}