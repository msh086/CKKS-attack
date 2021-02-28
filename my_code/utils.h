//
// Created by msh on 2021/1/23.
//

#ifndef CKKS_UTILS_H
#define CKKS_UTILS_H
#include <vector>
#include <complex>
#include <stdexcept>
#include <map>
#include <NTL/ZZ.h>

namespace utils {
    template<typename T>
    double max_error(std::vector<T>& baseline, std::vector<T>& noisy){
        int len = baseline.size();
        if(len != noisy.size())
            throw std::invalid_argument("array lengths do not match");
        double max_err = 0;
        for(int i = 0; i < len; i++){
            max_err = std::max(max_err, norm(baseline[i] - noisy[i]));
        }
        return max_err;
    }

    template <typename T>
    double relative_error(std::vector<T>& baseline, std::vector<T>& noisy){
        int len = baseline.size();
        if(len != noisy.size())
            throw std::invalid_argument("array lengths do not match");
        double max_rel_err = 0;
        for(int i = 0; i < len; i++){
            double norm_base = norm(baseline[i]);
            if(norm_base == 0)
                continue;
            double err = norm(baseline[i] - noisy[i]);
            max_rel_err = std::max(max_rel_err, err / norm_base);
        }
        return max_rel_err;
    }

    template<typename T>
    void element_wise_mult(const std::vector<T>& op1, const std::vector<T>& op2, std::vector<T>& dst){
        int len = op1.size();
        if(len != op2.size())
            throw std::invalid_argument("array lengths do not match");
        dst.resize(len);
        for(int i = 0; i < len; i++){
            dst[i] = op1[i] * op2[i];
        }
    }

    template<typename T>
    void element_wise_mult_inplace(std::vector<T>& multiplicand, const std::vector<T>& multiplier){
        element_wise_mult(multiplicand, multiplier, multiplicand);
    }

    template<typename T>
    void element_wise_add(const std::vector<T>& op1, const std::vector<T>& op2, std::vector<T>& dst){
        int len = op1.size();
        if(len != op2.size())
            throw std::invalid_argument("array lengths do not match");
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = op1[i] + op2[i];
    }

    template<typename T>
    void element_wise_add_inplace(std::vector<T>& addand, const std::vector<T>& adder){
        element_wise_add(addand, adder, addand);
    }

    template<typename T>
    void element_wise_substract(const std::vector<T>& op1, const std::vector<T>& op2, std::vector<T>& dst){
        int len = op1.size();
        if(len != op2.size())
            throw std::invalid_argument("array lengths do not match");
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = op1[i] - op2[i];
    }

    template<typename T>
    void element_wise_substract_inplace(std::vector<T>& suband, const std::vector<T>& subber){
        element_wise_substract(suband, subber, suband);
    }

    template<typename T>
    void element_wise_pow(const std::vector<T>& src, std::vector<T>& dst, int power){
        int len = src.size();
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = std::pow(src[i], power);
    }

    template<typename T>
    void element_wise_pow_inplace(std::vector<T>& src, int power){
        element_wise_pow(src, src, power);
    }

    template<typename T>
    void scale_mult(const std::vector<T>& src, std::vector<T>& dst, T factor){
        int len = src.size();
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = src[i] * factor;
    }

    template<typename T>
    void scale_mult_inplace(std::vector<T>& src, T factor){
        scale_mult(src, src, factor);
    }

    template<typename T>
    void scale_add(const std::vector<T>& src, std::vector<T>& dst, T factor){
        int len = src.size();
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = src[i] + factor;
    }

    template<typename T>
    void scale_add_inplace(std::vector<T>& src, T factor){
        scale_add(src, src, factor);
    }

    template<typename T>
    void to_norm(const std::vector<T>& src, std::vector<double>& dst){
        int len = src.size();
        dst.resize(len);
        for(int i = 0; i < len; i++)
            dst[i] = norm(src[i]);
    }

    template<typename T>
    double norm(std::complex<T> val){
        return std::sqrt(double(val.imag() * val.imag() + val.real() * val.real()));
    }

    template<typename T>
    double norm(T val){
        return std::abs(val);
    }

    template<typename T>
    T max_ele(const std::vector<T>& src){
        T res = 0;
        for(auto ele : src)
            res = std::max(ele, res);
        return res;
    }

    template<typename T>
    double infty_norm(std::vector<T>& vec) {
        double res = 0;
        for(auto val : vec)
            res = std::max(res, norm(val));
        return res;
    }

    double mean(std::vector<double>& vec);

    void sample_complex_circle(std::vector<std::complex<double>>& dst, size_t n, double radius);

    void sample_double(std::vector<double>& dst, size_t n, double radius);

    template<typename T>
    void element_wise_eval_polynomial(const std::vector<double>& coeffs, const std::vector<T>& src,
                                      std::vector<T>& dst){
        int len = src.size();
        dst.resize(len);
        std::fill_n(dst.begin(), len, 0);
        int n_coeffs = coeffs.size();
        if(n_coeffs == 0){
            printf("no plaintext polynomial evaluated\n");
            dst = src;
            return;
        }

        scale_add_inplace(dst, static_cast<T>(coeffs[0]));

        std::vector<T> pow_vec(src), tmp_vec;
        for(int i = 1; i < n_coeffs; i++){
            scale_mult(pow_vec, tmp_vec, static_cast<T>(coeffs[i]));
            element_wise_add_inplace(dst, tmp_vec);

            element_wise_mult_inplace(pow_vec, src);
        }
    }

    // note that `element_wise_eval_polynomial(coeffs, src, src)` is wrong
    template<typename T>
    void element_wise_eval_polynomial_inplace(const std::vector<double>& coeffs, std::vector<T>& src){
        std::vector<T> copy = src;
        element_wise_eval_polynomial(coeffs, copy, src);
    }

    // template classes & functions below are for fast calculation of total bits (or floor(log2(x)))

    template<typename T, unsigned int byte_size = sizeof(T)>
    struct valid_bits_helper{}; // never reach here

    template<typename T>
    struct valid_bits_helper<T, 8>{
        static constexpr unsigned int get(T val){
            return 64 - __builtin_clzl(val);
        }
    };

    template<typename T>
    struct valid_bits_helper<T, 4>{
        static constexpr unsigned int get(T val){
            return 32 - __builtin_clz(val);
        }
    };

    template<typename T>
    constexpr unsigned int count_valid_bits(T val){
        return valid_bits_helper<T>::get(val);
    }
    // end

    enum FuncType{
        F_NONE,
        F_SIGMOID,
        F_EXPONENT
    };

    extern const std::map<FuncType, std::vector<double>> func_map;

    // deprecated
    /*
    template<FuncType funcType>
    class MaclaurinSeries{
    public:
        static std::vector<double> gen_coeffs(int deg){
            return std::vector<double>{0., 1.};
        }
    };

    template<>
    class MaclaurinSeries<FuncType::F_EXPONENT>{
    public:
        static std::vector<double> gen_coeffs(int deg){
            std::vector<double> coeffs(deg + 1);
            coeffs[0] = 1;
            for(int i = 1; i <= deg; i++)
                coeffs[i] = coeffs[i] / i;
            return coeffs;
        }
    };

    template<>
    class MaclaurinSeries<FuncType::F_SIGMOID>{
        static const std::vector<double> precalced_coeffs;
    public:
        static std::vector<double> gen_coeffs(int deg){
            if(precalced_coeffs.size() < deg + 1)
                throw std::invalid_argument("Required degree greater than pre-calculated degree");
            std::vector<double> coeffs = precalced_coeffs;
            coeffs.resize(deg + 1);
            return coeffs;
        }
    };

    // reference to https://mathworld.wolfram.com/SigmoidFunction.html
    const std::vector<double> MaclaurinSeries<FuncType::F_SIGMOID>::precalced_coeffs = {
            1./2, 1./4, 0, -1./48, 0, 1./480, 0, -17./80640, 0, 31./1451520, 0
    };
    */
};


#endif //CKKS_UTILS_H
