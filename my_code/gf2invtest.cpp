//
// Created by msh on 2021/3/19.
// test the invertibility of R_2 elements
//

#include <NTL/GF2E.h>
#include <iostream>

constexpr uint64_t mask1 = 0x5555555555555555;
constexpr uint64_t mask2 = 0x3333333333333333;
constexpr uint64_t mask4 = 0x0f0f0f0f0f0f0f0f;
constexpr uint64_t mask8 = 0x00ff00ff00ff00ff;
constexpr uint64_t mask16 = 0x0000ffff0000ffff;
constexpr uint64_t mask32 = 0x00000000ffffffff;

uint64_t countOnes(uint64_t val) {
    val = (val & mask1) + ((val >> 1) & mask1); // 32 * 2 bits
    val = (val & mask2) + ((val >> 2) & mask2); // 16 * 4 bits
    val = (val & mask4) + ((val >> 4) & mask4); // 8 * 8 bits
    val = (val & mask8) + ((val >> 8) & mask8); // 4 * 16 bits
    val = (val & mask16) + ((val >> 16) & mask16); // 2 * 32 bits
    val = (val & mask32) + ((val >> 32) & mask32); // 1 * 64 bits
    return val;
}

uint64_t oddOnes(uint64_t val) {
    val ^= val >> 1;
    val ^= val >> 2;
    val ^= val >> 4;
    val ^= val >> 8;
    val ^= val >> 16;
    val ^= val >> 32;
    return val & 1;
}

int main() {
    uint64_t logn = 4;
    NTL::GF2X untruncated;
    uint64_t n = 1 << logn, cur = 1;
    uint64_t floor = uint64_t(1) << n; // Note: n < 64 or logn < 6 is required
    /**
     * logn = 3, 4, 5
     * n = 8, 16, 32
     * floor = 256, 65536, 2^32
     * */

    untruncated.SetLength(n + 1);
    untruncated[0] = 1;
    untruncated[n] = 1;
    NTL::GF2E::init(untruncated);

    int non_invertible_count = 0;
    while(cur < floor) {
        if(oddOnes(cur)) {
            untruncated.SetLength(n);
            for(int i = 0; i < n; i++) {
                if(cur & (1 << i))
                    untruncated[i] = 1;
                else
                    untruncated[i] = 0;
            }
            untruncated.normalize();
            try {
                NTL::inv(NTL::conv<NTL::GF2E>(untruncated));
            } catch (NTL::ErrorObject& e) {
                non_invertible_count++;
            }
        }
        cur++;
    }
    printf("there are %d non-invertible elements in R_2 that estimates to 1 at 1. log(n) = %lu\n", non_invertible_count, logn);
}