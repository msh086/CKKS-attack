//
// Created by msh on 2021/1/21.
//

#ifndef CKKS_SEALATTACK_H
#define CKKS_SEALATTACK_H
#include <seal/util/iterator.h>
#include <seal/context.h>
#include <vector>
#include <complex>

using namespace seal;

class sealattack {
public:
    static void add_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    static void sub_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    static void mult_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                  const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    static bool inv_ntt(util::ConstCoeffIter src, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    static void multi_precision_to_crt(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle& pool);

    static void crt_to_multi_precision(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle& pool);

    static void crt_to_ntt(util::CoeffIter src, const std::vector<Modulus>& coeff_modulus,
                           size_t coeff_count, const util::NTTTables* ntt_tables);

    static void ntt_to_crt(util::CoeffIter src, const std::vector<Modulus>& coeff_modulus,
                           size_t coeff_count, const util::NTTTables* ntt_tables);

    static void multi_precision_to_double_unscaled(util::ConstCoeffIter src, std::vector<double>& dst,
                                                   const SEALContext::ContextData &ctxt);

    static void sample_gaussian_noise(util::CoeffIter dst, const SEALContext::ContextData &ctxt, double std_dev, double max_dev);

};


#endif //CKKS_SEALATTACK_H
