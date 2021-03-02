//
// Created by msh on 2021/1/21.
//

#ifndef CKKS_SEALATTACK_H
#define CKKS_SEALATTACK_H
#include <seal/util/iterator.h>
#include <seal/context.h>
#include <vector>
#include <complex>
#include <seal/seal.h>

using namespace seal;

namespace sealattack {
    void add_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    void sub_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    void mult_ntt(util::ConstCoeffIter op1, util::ConstCoeffIter op2, util::CoeffIter dst,
                  const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    bool inv_ntt(util::ConstCoeffIter src, util::CoeffIter dst,
                 const std::vector<Modulus>& coeff_modulus, size_t coeff_count);

    void multi_precision_to_crt(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle& pool);

    void crt_to_multi_precision(util::CoeffIter src, const SEALContext::ContextData &ctxt, size_t coeff_count,
                                       MemoryPoolHandle& pool);

    void crt_to_ntt(util::CoeffIter src, const std::vector<Modulus>& coeff_modulus,
                           size_t coeff_count, const util::NTTTables* ntt_tables);

    void ntt_to_crt(util::CoeffIter src, const std::vector<Modulus>& coeff_modulus,
                           size_t coeff_count, const util::NTTTables* ntt_tables);

    void multi_precision_to_double_unscaled(util::ConstCoeffIter src, std::vector<double>& dst,
                                                   const SEALContext::ContextData &ctxt);

    void sample_gaussian_noise(util::CoeffIter dst, const SEALContext::ContextData &ctxt, double std_dev, double max_dev);


    /**
     * @brief evaluate a polynomial with floating point coefficients homomorphically
     *
     *
     * @param[in] coeffs: vector of coefficients, starting from low degree terms
     * @param[in] evaluator: evaluator that performs homomorphic operations
     * @param[in] relinKeys: relinearization keys
     * @param[in] x: input ciphertext, the variable in the to-be-evaluated polynomial
     * @param[out] dst: destination to store the result
     * */
    void evaluate_poly(const std::vector<double>& coeffs, seal::Evaluator& evaluator, seal::CKKSEncoder& encoder,
                       const seal::RelinKeys& relinKeys, const seal::Ciphertext& x, seal::Ciphertext& dst);

    void evaluate_poly_inplace(const std::vector<double>& coeffs, seal::Evaluator& evaluator, seal::CKKSEncoder& encoder,
                               const seal::RelinKeys& relinKeys, seal::Ciphertext& x);

    void calc_mean(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                   const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                   const seal::Ciphertext& src, seal::Ciphertext& dst);

    void calc_mean_inplace(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                           const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                           seal::Ciphertext& src);

    void calc_variance(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                       const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                       const seal::Ciphertext& src, seal::Ciphertext& dst, const seal::Ciphertext *mean = nullptr);

    void calc_variance_inplace(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                               const seal::GaloisKeys& galoisKeys, const seal::RelinKeys& relinKeys,
                               seal::Ciphertext& src, const seal::Ciphertext *mean = nullptr);

    void scale_to_same_level(seal::Evaluator& evaluator, seal::CKKSEncoder& encoder, const seal::SEALContext& ctxt,
                             const seal::RelinKeys& relinKeys, seal::Ciphertext& c1, seal::Ciphertext& c2);
};


#endif //CKKS_SEALATTACK_H
