//
// Created by msh on 2021/1/21.
//

#include "sealattack.h"
#include <seal/util/polyarithsmallmod.h>
#include <seal/randomtostd.h>
#include <seal/randomgen.h>
#include <seal/util/clipnormal.h>
#include <seal/util/iterator.h>
#include "utils.h"

namespace sealattack {
    void add_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                             const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
        auto coeff_modulus_count = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_count; i++) {
            util::add_poly_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count, coeff_modulus[i],
                                    dst + (i * coeff_count));
        }
    }

    void mult_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                              const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
        auto coeff_modulus_count = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_count; i++) {
            util::dyadic_product_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count,
                                          coeff_modulus[i],
                                          dst + (i * coeff_count));
        }
    }

    void sub_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                             const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
        auto coeff_modulus_count = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_count; i++) {
            util::sub_poly_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count, coeff_modulus[i],
                                    dst + (i * coeff_count));
        }
    }

    bool inv_ntt(util::ConstCoeffIter src, util::CoeffIter dst, const std::vector<Modulus> &coeff_modulus,
                             size_t coeff_count) {
        auto coeff_modulus_count = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_count; i++) {
            for (size_t j = 0; j < coeff_count; j++) {
                if (!util::try_invert_uint_mod(*(src + (i * coeff_count + j)), coeff_modulus[i],
                                               *(dst + (i * coeff_count + j))))
                    return false;
            }
        }
        return true;
    }

    void crt_to_ntt(util::CoeffIter src, const std::vector<Modulus> &coeff_modulus,
                                size_t coeff_count, const util::NTTTables *ntt_tables) {
        auto coeff_modulus_size = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_size; i++) {
            util::ntt_negacyclic_harvey(src + (i * coeff_count), ntt_tables[i]);
        }
    }

    void ntt_to_crt(util::CoeffIter src, const std::vector<Modulus> &coeff_modulus,
                                size_t coeff_count, const util::NTTTables *ntt_tables) {
        auto coeff_modulus_size = coeff_modulus.size();
        for (size_t i = 0; i < coeff_modulus_size; i++) {
            util::inverse_ntt_negacyclic_harvey(src + (i * coeff_count), ntt_tables[i]);
        }
    }

    void
    multi_precision_to_crt(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle &pool) {
        ctxt.rns_tool()->base_q()->decompose_array(src, coeff_count, pool);
    }

    void
    crt_to_multi_precision(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle &pool) {
        ctxt.rns_tool()->base_q()->compose_array(src, coeff_count, pool);
    }

    void multi_precision_to_double_unscaled(util::ConstCoeffIter src, std::vector<double> &dst,
                                                        const SEALContext::ContextData &ctxt) {
        auto modulus_product = ctxt.total_coeff_modulus();
        auto upper_half_threshold = ctxt.upper_half_threshold();
        auto &params = ctxt.parms();
        auto coeff_count = params.poly_modulus_degree();
        auto coeff_modulus_size = params.coeff_modulus().size();

        // Create floating-point representations of the multi-precision integer coefficients
        double two_pow_64 = std::pow(2.0, 64);
        dst.resize(coeff_count);
        for (std::size_t i = 0; i < coeff_count; i++) {
            dst[i] = 0.0;
            if (util::is_greater_than_or_equal_uint(
                    src + (i * coeff_modulus_size), upper_half_threshold, coeff_modulus_size)) {
                double scaled_two_pow_64 = 1.0;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64) {
                    if (src[i * coeff_modulus_size + j] > modulus_product[j]) {
                        auto diff = src[i * coeff_modulus_size + j] - modulus_product[j];
                        dst[i] += diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    } else {
                        auto diff = modulus_product[j] - src[i * coeff_modulus_size + j];
                        dst[i] -= diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                    }
                }
            } else {
                double scaled_two_pow_64 = 1.0;
                for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64) {
                    auto curr_coeff = src[i * coeff_modulus_size + j];
                    dst[i] += curr_coeff ? static_cast<double>(curr_coeff) * scaled_two_pow_64 : 0.0;
                }
            }

            // Scaling instead incorporated above; this can help in cases
            // where otherwise pow(two_pow_64, j) would overflow due to very
            // large coeff_modulus_size and very large scale
            // dst[i] = res_accum * inv_scale;
        }
    }

    void sample_gaussian_noise(util::CoeffIter dst, const SEALContext::ContextData &ctxt, double std_dev,
                                           double max_dev) {
        auto &params = ctxt.parms();
        auto prng = params.random_generator()->create();
        auto &coeff_modulus = params.coeff_modulus();
        size_t coeff_modulus_size = coeff_modulus.size();
        size_t coeff_count = params.poly_modulus_degree();
        auto ntt_tables = ctxt.small_ntt_tables();

        RandomToStandardAdapter engine(prng);
        util::ClippedNormalDistribution dist(0, std_dev, max_dev);
        SEAL_ITERATE(dst, coeff_count, [&](auto &I) {
            auto noise = static_cast<int64_t>(dist(engine));
            auto flag = static_cast<uint64_t>(-static_cast<int64_t>(noise < 0));
            SEAL_ITERATE(
                    iter(util::StrideIter<uint64_t *>(&I, coeff_count), coeff_modulus), coeff_modulus_size,
                    [&](auto J) { *std::get<0>(J) = static_cast<uint64_t>(noise) + (flag & std::get<1>(J).value()); });
        });
        crt_to_ntt(dst, coeff_modulus, coeff_count, ntt_tables);
    }

    void evaluate_poly(const std::vector<double>& coeffs, seal::Evaluator& evaluator, seal::CKKSEncoder& encoder,
                       const seal::RelinKeys& relinKeys, const seal::Ciphertext& x, seal::Ciphertext& dst){
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
        std::vector<seal::Ciphertext> tower(max_height + 1);
        tower[0] = x;
        for(int i = 1; i <= max_height; i++){
            evaluator.square(tower[i - 1], tower[i]);
            evaluator.relinearize_inplace(tower[i], relinKeys);
            evaluator.rescale_to_next_inplace(tower[i]);
        }

        // initialize dst, since we are not allowed to set dst to 0 and add all terms up
        seal::Plaintext tmp_plaintext;
        dst = x;
        encoder.encode(coeffs[1], dst.parms_id(), dst.scale(), tmp_plaintext);
        evaluator.multiply_plain_inplace(dst, tmp_plaintext);
        evaluator.relinearize_inplace(dst, relinKeys);
        evaluator.rescale_to_next_inplace(dst);
        // now dst = c1 * x
        encoder.encode(coeffs[0], dst.parms_id(), dst.scale(), tmp_plaintext);
        evaluator.add_plain_inplace(dst, tmp_plaintext);
        // now dst = c0 + c1 * x

        // a tmp_plaintext encoding a vector of ones enables us to rescale ciphertext at higher levels to lower levels
        // an alternative is to change the scales to the same value, then mod switch to the lower level
        std::vector<seal::Plaintext> plaintexts_all_ones(max_height + 1);
        int ciphertext_level = 1; // 1 rescaling is taken
        for(int i = 0; i <= max_height; i++)
            encoder.encode(1, tower[i].parms_id(), tower[i].scale(), plaintexts_all_ones[i]);
        // then sum up other terms
        seal::Ciphertext tmp_ciphertext; // for scratch
        for(int cur_deg = 2; cur_deg < n_coeffs; cur_deg++){

            unsigned int cur_deg_total_bits = utils::count_valid_bits(cur_deg), cursor_bit_idx = 0;
            for(; cursor_bit_idx < cur_deg_total_bits; cursor_bit_idx++) {
                if ((1 << cursor_bit_idx) & cur_deg)
                    break;
            }

            if (fabs(coeffs[cur_deg]) * tower[cursor_bit_idx].scale() < 0.5) // too small s.t. encoding results is zero poly
                continue;

            tmp_ciphertext = tower[cursor_bit_idx];
            encoder.encode(coeffs[cur_deg], tmp_ciphertext.parms_id(), tmp_ciphertext.scale(), tmp_plaintext);
            evaluator.multiply_plain_inplace(tmp_ciphertext, tmp_plaintext);
            evaluator.relinearize_inplace(tmp_ciphertext, relinKeys);
            evaluator.rescale_to_next_inplace(tmp_ciphertext);
            while(++cursor_bit_idx < cur_deg_total_bits){
                if((1 << cursor_bit_idx) & cur_deg){
                    evaluator.multiply_inplace(tmp_ciphertext, tower[cursor_bit_idx]);
                    evaluator.relinearize_inplace(tmp_ciphertext, relinKeys);
                    evaluator.rescale_to_next_inplace(tmp_ciphertext);
                }
                else{
                    evaluator.multiply_plain_inplace(tmp_ciphertext, plaintexts_all_ones[cursor_bit_idx]);
                    evaluator.relinearize_inplace(tmp_ciphertext, relinKeys);
                    evaluator.rescale_to_next_inplace(tmp_ciphertext);
                }
            }
            // rescale dst to the same level as tmp_ciphertext, then add them together
            // level of tmp_ciphertext is guaranteed to be cur_deg_total_bits
            // the level of dst need to be recorded to perform rescaling
            while(dst.parms_id() != tmp_ciphertext.parms_id()) {
                evaluator.multiply_plain_inplace(dst, plaintexts_all_ones[ciphertext_level++]);
                evaluator.relinearize_inplace(dst, relinKeys);
                evaluator.rescale_to_next_inplace(dst);
            }
            evaluator.add_inplace(dst, tmp_ciphertext);
        }
        // the homomorphic calculation will raise the level of the ciphertext by count_valid_bits(n_coeffs - 1)
        // or bits(max_degree)
    }

    void evaluate_poly_inplace(const std::vector<double>& coeffs, seal::Evaluator& evaluator, seal::CKKSEncoder& encoder,
                               const seal::RelinKeys& relinKeys, seal::Ciphertext& x){
        seal::Ciphertext copy = x;
        evaluate_poly(coeffs, evaluator, encoder, relinKeys, copy, x);
    }
}