//
// Created by msh on 2021/1/21.
//

#include "sealattack.h"
#include <seal/util/polyarithsmallmod.h>
#include <seal/randomtostd.h>
#include <seal/randomgen.h>
#include <seal/util/clipnormal.h>
#include <seal/util/iterator.h>
void sealattack::add_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                         const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
    auto coeff_modulus_count = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_count; i++){
        util::add_poly_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count, coeff_modulus[i],
                                dst + (i * coeff_count));
    }
}

void sealattack::mult_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                          const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
    auto coeff_modulus_count = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_count; i++){
        util::dyadic_product_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count, coeff_modulus[i],
                                      dst + (i * coeff_count));
    }
}

void sealattack::sub_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                         const std::vector<Modulus> &coeff_modulus, size_t coeff_count) {
    auto coeff_modulus_count = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_count; i++){
        util::sub_poly_coeffmod(op1 + (i * coeff_count), op2 + (i * coeff_count), coeff_count, coeff_modulus[i],
                                      dst + (i * coeff_count));
    }
}

bool sealattack::inv_ntt(util::ConstCoeffIter src, util::CoeffIter dst, const std::vector<Modulus> &coeff_modulus,
                         size_t coeff_count) {
    auto coeff_modulus_count = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_count; i++){
        for(size_t j = 0; j < coeff_count; j++){
            if(!util::try_invert_uint_mod(*(src + (i * coeff_count + j)), coeff_modulus[i], *(dst + (i * coeff_count + j))))
                return false;
        }
    }
    return true;
}

void sealattack::crt_to_ntt(util::CoeffIter src, const std::vector<Modulus> &coeff_modulus,
                            size_t coeff_count, const util::NTTTables *ntt_tables) {
    auto coeff_modulus_size = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_size; i++){
        util::ntt_negacyclic_harvey(src + (i * coeff_count), ntt_tables[i]);
    }
}

void sealattack::ntt_to_crt(util::CoeffIter src, const std::vector<Modulus> &coeff_modulus,
                            size_t coeff_count, const util::NTTTables *ntt_tables) {
    auto coeff_modulus_size = coeff_modulus.size();
    for(size_t i = 0; i < coeff_modulus_size; i++){
        util::inverse_ntt_negacyclic_harvey(src + (i * coeff_count), ntt_tables[i]);
    }
}

void sealattack::multi_precision_to_crt(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                        MemoryPoolHandle &pool) {
    ctxt.rns_tool()->base_q()->decompose_array(src, coeff_count, pool);
}

void sealattack::crt_to_multi_precision(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                        MemoryPoolHandle& pool) {
    ctxt.rns_tool()->base_q()->compose_array(src, coeff_count, pool);
}

void sealattack::multi_precision_to_double_unscaled(util::ConstCoeffIter src, std::vector<double>& dst,
                                                    const SEALContext::ContextData &ctxt) {
    auto modulus_product = ctxt.total_coeff_modulus();
    auto upper_half_threshold = ctxt.upper_half_threshold();
    auto &params = ctxt.parms();
    auto coeff_count = params.poly_modulus_degree();
    auto coeff_modulus_size = params.coeff_modulus().size();

    // Create floating-point representations of the multi-precision integer coefficients
    double two_pow_64 = std::pow(2.0, 64);
    dst.resize(coeff_count);
    for (std::size_t i = 0; i < coeff_count; i++)
    {
        dst[i] = 0.0;
        if (util::is_greater_than_or_equal_uint(
                src + (i * coeff_modulus_size), upper_half_threshold, coeff_modulus_size))
        {
            double scaled_two_pow_64 = 1.0;
            for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
            {
                if (src[i * coeff_modulus_size + j] > modulus_product[j])
                {
                    auto diff = src[i * coeff_modulus_size + j] - modulus_product[j];
                    dst[i] += diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                }
                else
                {
                    auto diff = modulus_product[j] - src[i * coeff_modulus_size + j];
                    dst[i] -= diff ? static_cast<double>(diff) * scaled_two_pow_64 : 0.0;
                }
            }
        }
        else
        {
            double scaled_two_pow_64 = 1.0;
            for (std::size_t j = 0; j < coeff_modulus_size; j++, scaled_two_pow_64 *= two_pow_64)
            {
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

void sealattack::sample_gaussian_noise(util::CoeffIter dst, const SEALContext::ContextData &ctxt, double std_dev,
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