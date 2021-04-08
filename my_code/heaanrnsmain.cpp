//
// Created by msh on 2021/4/8.
//

#include "TestScheme.h"
#include "Numb.h"
#include "Context.h"
#include "SecretKey.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "SchemeAlgo.h"

#include "utils.h"
#include <vector>
#include <complex>

int main(){
    long logN = 16, logp = 40, levels = 5, specials = 1;
    long slots = 1 << (logN - 1);
    Context context(logN, logp, levels, specials);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    using cmpl64 = std::complex<double>;
    std::vector<cmpl64> plain_vec;
    utils::sample_complex_circle(plain_vec, slots, 1);

    Ciphertext cipher = scheme.encrypt(plain_vec.data(), slots, levels);
    auto dec_res = scheme.decrypt(secretKey, cipher);

    std::vector<cmpl64> dec_vec(dec_res, dec_res + slots);
    delete[] dec_res;

    StringUtils::showcompare(plain_vec.data(), dec_vec.data(), 6, "val");

    return 0;
}