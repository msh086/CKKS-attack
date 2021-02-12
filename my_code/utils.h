//
// Created by msh on 2021/1/23.
//

#ifndef CKKS_UTILS_H
#define CKKS_UTILS_H
#include <vector>
#include <complex>
#include <stdexcept>
#include <initializer_list>
#include <NTL/ZZ.h>

class utils {
public:
    template<typename T>
    static double max_error(std::vector<T>& baseline, std::vector<T>& noisy){
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
    static double relative_error(std::vector<T>& baseline, std::vector<T>& noisy){
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
    static void element_wise_mult_inplace(std::vector<T>& multiplicand, const std::vector<T>& multiplier){
        element_wise_mult(multiplicand, multiplier, multiplicand);
    }

    template<typename T>
    static void element_wise_mult(const std::vector<T>& op1, const std::vector<T>& op2, std::vector<T>& dst){
        int len = op1.size();
        if(len != op2.size())
            throw std::invalid_argument("array lengths do not match");
        dst.resize(len);
        for(int i = 0; i < len; i++){
            dst[i] = op1[i] * op2[i];
        }
    }

    template<typename T>
    static double norm(std::complex<T> val){
        return std::sqrt(double(val.imag() * val.imag() + val.real() * val.real()));
    }

    template<typename T>
    static double norm(T val){
        return std::abs(val);
    }

    template<typename T>
    static double infty_norm(std::vector<T>& vec) {
        double res = 0;
        for(auto val : vec)
            res = std::max(res, norm(val));
        return res;
    }

    static double mean(std::vector<double>& vec);

    static void sample_complex_circle(std::vector<std::complex<double>>& dst, size_t n, double radius);

    static void sample_double(std::vector<double>& dst, size_t n, double radius);
};


#endif //CKKS_UTILS_H
