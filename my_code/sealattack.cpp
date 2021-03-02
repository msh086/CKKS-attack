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
        if(&x == &dst){
            throw std::invalid_argument("src and dst Ciphertext should not have the same address");
        }
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

    void calc_variance(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                       const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                       const seal::Ciphertext& src, seal::Ciphertext& dst, const seal::Ciphertext *mean){
        if(&src == &dst || mean == &dst)
            throw std::invalid_argument("dst cannot have the same address as src or mean");
        seal::Ciphertext tmp_ciphertext;
        evaluator.complex_conjugate(src, galoisKeys, tmp_ciphertext);
        evaluator.multiply_inplace(tmp_ciphertext, src);
        evaluator.relinearize_inplace(tmp_ciphertext, relinKeys);
        evaluator.rescale_to_next_inplace(tmp_ciphertext);
        // compute mean(x^2)
        calc_mean(evaluator, encoder, ctxt, galoisKeys, relinKeys, tmp_ciphertext, dst);
        // compute mean(x)^2
        if(mean)
            tmp_ciphertext = *mean;
        else
            calc_mean(evaluator, encoder, ctxt, galoisKeys, relinKeys, src, tmp_ciphertext);
        seal::Ciphertext tmp_ciphertext_2;
        evaluator.complex_conjugate(tmp_ciphertext, galoisKeys, tmp_ciphertext_2);
        evaluator.multiply_inplace(tmp_ciphertext, tmp_ciphertext_2);
        evaluator.relinearize_inplace(tmp_ciphertext, relinKeys);
        evaluator.rescale_to_next_inplace(tmp_ciphertext);
        // compute var(x) = mean(x^2) - mean(x)^2
        scale_to_same_level(evaluator, encoder, ctxt, relinKeys, dst, tmp_ciphertext);
        evaluator.sub_inplace(dst, tmp_ciphertext);
    }

    void calc_variance_inplace(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                               const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                               seal::Ciphertext& src, const seal::Ciphertext *mean){
        seal::Ciphertext src_copy = src;
        calc_variance(evaluator, encoder, ctxt, galoisKeys, relinKeys, src_copy, src, mean);
    }

    void calc_mean(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                   const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                   const seal::Ciphertext& src, seal::Ciphertext& dst){
        if(&src == &dst)
            throw std::invalid_argument("src and dst should not have the same address");
        dst = src;
        // n_slot is guaranteed to be power of 2
        int n_slot = ctxt.get_context_data(src.parms_id())->parms().poly_modulus_degree() >> 1;
        int n_slot_copy = n_slot;
        int step = 1;
        seal::Ciphertext rotated_tmp;
        while(n_slot_copy>>=1){
            evaluator.rotate_vector(dst, step, galoisKeys, rotated_tmp);
            evaluator.add_inplace(dst, rotated_tmp);
            step <<= 1;
        }
        seal::Plaintext plain_tmp;
        encoder.encode(1. / n_slot, dst.parms_id(), dst.scale(), plain_tmp);
        evaluator.multiply_plain_inplace(dst, plain_tmp);
        evaluator.relinearize_inplace(dst, relinKeys);
        evaluator.rescale_to_next_inplace(dst);
    }

    void calc_mean_inplace(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                   const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                   seal::Ciphertext& src){
        seal::Ciphertext src_copy = src;
        calc_mean(evaluator, encoder, ctxt, galoisKeys, relinKeys, src_copy, src);
    }

    void scale_to_same_level(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                             const seal::RelinKeys& relinKeys, seal::Ciphertext& c1, seal::Ciphertext& c2){
        auto remain1 = ctxt.get_context_data(c1.parms_id())->parms().coeff_modulus().size(),
            remain2 = ctxt.get_context_data(c2.parms_id())->parms().coeff_modulus().size();
        seal::Ciphertext* higher = nullptr;
        std::size_t diff = 0;
        if(remain1 < remain2){
            diff = remain2 - remain1;
            higher = &c2;
        }
        else{
            diff = remain2 - remain1;
            higher = &c1;
        }
        seal::Plaintext tmp_plain;
        while(diff){ // using `diff--` here may lead to underflow
            encoder.encode(1, higher->parms_id(), higher->scale(), tmp_plain);
            evaluator.multiply_plain_inplace(*higher, tmp_plain);
            evaluator.relinearize_inplace(*higher, relinKeys);
            evaluator.rescale_to_next_inplace(*higher);
            diff--;
        }
    }
}