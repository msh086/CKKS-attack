//
// Created by msh on 2021/3/2.
//

#include <memory>
#include <vector>
#include "utils.h"
#include "helibattack.h"
#include <NTL/ZZ.h>
#include <iostream>
namespace helibattack {
    void eval_polynomial(const std::vector<double>& coeffs, const helib::Ctxt& x, helib::Ctxt& dst){
        int n_coeffs = coeffs.size();
        // trivial cases
        if(n_coeffs == 0) {
            printf("no homomorphic evaluation\n");
            dst = x;
            return;
        }
        if(n_coeffs == 1){
            throw std::invalid_argument("homomorphic polynomial is constant\n");
        }

        unsigned int max_height = utils::count_valid_bits(n_coeffs - 1) - 1;
        // max height is bits(max_degree) - 1

        // now generate x, x^2, ..., x^(2^max_height)
        std::vector<helib::Ctxt> tower;
        tower.reserve(max_height + 1);
        tower.emplace_back(x);
        for(int i = 1; i <= max_height; i++){
            tower.emplace_back(tower[i - 1]);
            tower[i].square();
            // NOTE: drop special primes, the 'P' in relinearization
            tower[i].dropSmallAndSpecialPrimes();
        }

        dst = x;
        dst *= coeffs[1];
        dst += coeffs[0];
        // now dst = c0 + c1 * x
        // then sum up other terms
        helib::Ctxt tmp_ciphertext(dst.getPubKey()); // for scratch
        for(int cur_deg = 2; cur_deg < n_coeffs; cur_deg++){

            unsigned int cur_deg_total_bits = utils::count_valid_bits(cur_deg), cursor_bit_idx = 0;
            for(; cursor_bit_idx < cur_deg_total_bits; cursor_bit_idx++) {
                if ((1 << cursor_bit_idx) & cur_deg)
                    break;
            }

            if (fabs(coeffs[cur_deg]) * tower[cursor_bit_idx].getRatFactor() < 0.5) // too small s.t. encoding results is zero poly
                continue;

            tmp_ciphertext = tower[cursor_bit_idx];
            tmp_ciphertext *= coeffs[cur_deg];
            while(++cursor_bit_idx < cur_deg_total_bits){
                if((1 << cursor_bit_idx) & cur_deg){
                    tmp_ciphertext *= tower[cursor_bit_idx];
                    // NOTE: drop special primes, the 'P' in relinearization
                    tmp_ciphertext.dropSmallAndSpecialPrimes();
                }
            }
            // TODO: addition between ciphertexts under different modulus q1, q2
            //  will up-scale both into modulus lcm(q1, q2),
            //  then, if the two scaling factors are not the same,
            //  another process is applied to make them match
            dst += tmp_ciphertext;
        }
    }

    void eval_polynomial_inplace(const std::vector<double>& coeffs, helib::Ctxt& x){
        eval_polynomial(coeffs, x, x);
    }

    void print_NTL_macros(){
        std::cout << "NTL_XD_BOUND: " << NTL_XD_BOUND << std::endl
                  << "NTL_XD_HBOUND: " << NTL_XD_HBOUND << std::endl
                  << "NTL_DOUBLE_PRECISION: " << NTL_DOUBLE_PRECISION << std::endl
                  << "NTL_BITS_PER_LONG: " << NTL_BITS_PER_LONG << std::endl;
    }

    void print_xdouble(NTL::xdouble val){
        std::cout << val << std::endl;
    }
}

