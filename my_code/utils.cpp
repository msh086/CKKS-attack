//
// Created by msh on 2021/1/23.
//

#include "utils.h"

double utils::mean(std::vector<double> &vec){
    double res = 0;
    for(auto val : vec)
        res += val;
    return res / vec.size();
}

void utils::sample_complex_circle(std::vector<std::complex<double>>& dst, size_t n, double radius){
    if(radius <= 0)
        throw std::invalid_argument("radius must be positive");
    dst.resize(n);
    double two_times_pi = 2 * M_PI;
    double r_mult = radius / double(0x10000), theta_mult = two_times_pi / double(0x100000000);
    for(auto &ele : dst){
        unsigned long long_bits = NTL::RandomWord();
        double r = std::sqrt(double(long_bits & 0xffffffff)) * r_mult;
        double theta = double(long_bits >> 32) * theta_mult;
        ele = std::polar(r, theta);
    }
}

void utils::sample_double(std::vector<double>& dst, size_t n, double radius){
    if(radius <= 0)
        throw std::invalid_argument("radius must be positive");
    dst.resize(n);
    for(auto &ele : dst){
        unsigned long long_bits = NTL::RandomWord();
        ele = double(long_bits >> 1) / double(0x8000000000000000) * ((long_bits & 1) ? 1 : -1) * radius;
    }
}