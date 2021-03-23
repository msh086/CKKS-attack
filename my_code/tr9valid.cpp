//
// Created by msh on 2021/3/23.
// given modulus q = 2^a, dimension n = 2^b and ring R_q = Z_q[X]/(X^n + 1)
// R_q* (the group of the units in R_q) consists exclusively of (S + 2T)
// where S is an unit in R_2 and T is an arbitrary polynomial in R_(2^(a-1))
// the program tests if NTL can compute the inverse correctly
//

#include <iostream>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>

uint64_t oddOnes(uint64_t val) {
    val ^= val >> 1;
    val ^= val >> 2;
    val ^= val >> 4;
    val ^= val >> 8;
    val ^= val >> 16;
    val ^= val >> 32;
    return val & 1;
}

constexpr uint64_t lowerBitsMask(uint32_t width) {
    return (uint64_t(1) << width) - 1;
}

template<uint32_t width, typename = typename std::enable_if<(width > 0 && width <= 63)>::type>
NTL::ZZX bits2ZZX(uint64_t bits, uint32_t maxDeg) {
    auto mask = lowerBitsMask(width);
    NTL::ZZX poly;
    poly.SetLength(maxDeg + 1);
    uint32_t coeff_idx = 0;
    while (bits && coeff_idx <= maxDeg) {
        poly[coeff_idx++] = bits & mask;
        bits >>= width;
    }
    poly.normalize();
    return poly;
}

bool invertible(const NTL::ZZX &op) {
    auto copy = op;
    copy.normalize();
    try {
        auto res = NTL::inv(NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(copy)));
    } catch (NTL::ErrorObject &e) {
        return false;
    }
    return true;
}

int main() {
    // total inverse count in R_2: 2^(n - 1)
    // total element in R_(2^(a - 1)): 2^((a - 1) * n)
    // total inverse test count: 2^(n - 1) * 2^((a - 1) * n) = 2^(n * a - 1)
    // NOTE: don't let n * a - 1 surpass 30, or it may be quite slow
    // total bits needed to store a polynomial in R_q is (a * n)
    // total bits needed to store a polynomial in R_2 is (n)
    constexpr uint64_t log2_q = 5, log2_n = 2;
    uint64_t q = uint64_t(1) << log2_q, n = uint64_t(1) << log2_n;
    uint64_t totalBits = ((uint64_t) 1 << log2_n) * log2_q;
    if (totalBits >= 63)
        return 1;
    // init R_q
    NTL::ZZ_p::init(NTL::ZZ(q));
    {
        NTL::ZZX modulus;
        modulus.SetLength(n + 1);
        modulus[0] = modulus[n] = 1;
        NTL::ZZ_pE::init(NTL::conv<NTL::ZZ_pX>(modulus));
    }
    // precompute all units in R_2
    uint64_t n_r2_units = 1 << (n - 1), r2_card = 1 << n;
    auto *r2_units = new NTL::ZZX[n_r2_units];
    for (uint64_t r2_poly = 0, cur = 0; r2_poly < r2_card; r2_poly++) {
        if (oddOnes(r2_poly)) {
            if (cur >= n_r2_units)
                return 1;
            r2_units[cur++] = bits2ZZX<1>(r2_poly, n - 1);
        }
    }
    // iterate through all R_(q/2) elements
    uint64_t rhq_poly_ceil = uint64_t(1) << (n * (log2_q - 1));
    for (uint64_t rhq_poly = 0; rhq_poly < rhq_poly_ceil; rhq_poly++) {
        auto rhq_element = bits2ZZX<log2_q>(rhq_poly, n - 1);
        // check all units in R_2
        for (int i = 0; i < n_r2_units; i++) {
            if (!invertible(2 * rhq_element + r2_units[i])) {
                std::cout << rhq_element << std::endl << r2_units[i] << std::endl;
//                return 1;
            }
        }
    }
    return 0;
}