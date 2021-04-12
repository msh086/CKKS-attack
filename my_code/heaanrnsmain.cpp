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
#include <NTL/ZZ.h>
#include <iostream>
#include "tclap/CmdLine.h"

using cmpl64 = std::complex<double>;

bool ntt_inverse(Context &ctxt, uint64_t *dst, const uint64_t *src, long levels) {
    const uint64_t *crt_part = src;
    for (int i = 0; i < levels; i++, crt_part += ctxt.N, dst += ctxt.N) {
        for (int j = 0; j < ctxt.N; j++) {
            if (crt_part[j] == 0)
                return false;
            // implicit type conversion from uint64_t to long, watch out for possible overflow
            dst[j] = NTL::InvMod(crt_part[j], ctxt.qVec[i]);
            // this never throws
        }
    }
    return true;
}

int modulus_bits(Context &ctxt) {
    int res = 0;
    for (int i = 0; i < ctxt.L; i++)
        res += NTL::NumBits(ctxt.qVec[i]);
    return res;
}

Plaintext my_decrypt_msg(Context &ctxt, const SecretKey &secretKey, const Ciphertext &ciphertext) {
    auto *plain = new uint64_t[ciphertext.l << ctxt.logN]();
    ctxt.mul(plain, ciphertext.ax, secretKey.sx, ciphertext.l);
    ctxt.add(plain, plain, ciphertext.bx, ciphertext.l);
    return Plaintext(plain, ciphertext.N, ciphertext.slots, ciphertext.l);
}

/**
 * memory layout:
 * [N] [N] ... [N], number of [N] is l
 * */
double *crt_compose(Context &ctxt, const uint64_t *src, long l) {
    // preprocessing
    NTL::ZZ productAll(1); // product of all primes
    std::vector<NTL::ZZ> primes(l);
    // arr[i] = P_i * (inverse of P_i wrt primes[i]), where P_i = productAll / primes[i]
    std::vector<NTL::ZZ> helperArr(l);
    NTL::ZZ tmp;
    for (int i = 0; i < l; i++) {
        primes[i] = ctxt.qVec[i];
        productAll *= primes[i];
    }
    NTL::ZZ halfModulus = productAll;
    halfModulus >>= 1;
    for (int i = 0; i < l; i++) {
        NTL::div(tmp, productAll, primes[i]); // tmp = P_i
        helperArr[i] = tmp % primes[i]; // arr[i] = P_i % primes[i]
        NTL::InvMod(helperArr[i], helperArr[i], primes[i]); // arr[i] = P_i^-1
        helperArr[i] *= tmp; // (P_i^-1 mod primes[i]) * P_i
    }
    // composing
    auto *res = new double[ctxt.N];
    const uint64_t *series;
    long maxBits = 0;
    for (int i = 0; i < ctxt.N; i++) {
        series = src + i;
        tmp = 0;
        for (int j = 0; j < l; j++, series += ctxt.N)
            NTL::AddMod(tmp, tmp, helperArr[j] * (*series), productAll);
        tmp %= productAll;
        if (tmp >= halfModulus)
            tmp -= productAll;
        // FIXME debug
        if (tmp % primes[0] != src[i])
            std::cout << "?";
        // end FIXME
        maxBits = std::max(maxBits, NTL::NumBits(tmp));
        res[i] = NTL::conv<double>(tmp);
    }
    std::cout << "max bits in coeffs: " << maxBits << std::endl;
    return res;
}


/**
 * copied from Context.cpp Context::decode line 503
 * the original version assumes that the coefficients fits in the first prime(about 60 bits)
 * here we use all the CRT primes to reconstruct a multi-precision integer(in crt_compose)
 * this is not necessary if the coeffs has significant bits less than 62
 *
 * NOTE: besides, the original version contains the conversion of the first prime into 'double' type
 *  this will introduce an error of magnitude exactly
 *  (consider the bit pattern of the first prime:
 *  1x...x0...01, with 62 bits in total and k-1 successive 0-bits ahead of the lowest 1-bit
 *  the highest 53 bits are kept, and the lower 9 bits are discarded
 *  which have the pattern 0...01(since k-1 is typically greater than 8)
 *  BUT, the coeff is also converted to double, then subtracted by the first prime(this conversion itself is accurate)
 *  following the rules for floating point arithmetic, the operand with a smaller exponent is scaled up so that its
 *  exponent matches the other,
 *  during this process, the lower 9 bits of the integral part are discarded,
 *  and a literally huge error is introduced, since those bits can take any value
 * */
void my_decode(Context &ctxt, const uint64_t *a, complex<double> *v, long slots, long l) {
    uint64_t *tmp = new uint64_t[l << ctxt.logN]();
    copy(a, a + (l << ctxt.logN), tmp);
    long gap = ctxt.Nh / slots;
    for (int i = 0; i < l; i++)
        ctxt.qiINTTAndEqual(tmp + (i << ctxt.logN), i);
    // reconstruct from CRT
    double *composed = crt_compose(ctxt, tmp, l);

    for (long j = 0, jdx = ctxt.Nh, idx = 0; j < slots; ++j, jdx += gap, idx += gap) {
        double mir = composed[idx] / ctxt.p;
        double mii = composed[jdx] / ctxt.p;
        v[j].real(mir);
        v[j].imag(mii);
    }
    ctxt.fftSpecial(v, slots);
}

cmpl64 *my_decode_alloc(Context &context, Plaintext &plain) {
    auto *res = new cmpl64[plain.slots];
    my_decode(context, plain.mx, res, plain.slots, plain.l);
    return res;
}

/**
 * the original encoding performs truncating instead of rounding, which introduces extra error
 * in this case, the re-encoded plaintext is identical as the original one iff.
 * the encoding error lies in [0, 1)
 *
 * I changed the behaviour to rounding here
 * */
void my_encode(Context &ctxt, uint64_t *a, const complex<double> *v, long slots, long l) {
    auto *uvals = new complex<double>[slots]();
    copy(v, v + slots, uvals);

    long gap = ctxt.Nh / slots;

    ctxt.fftSpecialInv(uvals, slots);

    for (long j = 0; j < slots; ++j) {
        uvals[j] *= ctxt.p;
    }

    for (long i = 0; i < l; ++i) {
        uint64_t *mi = a + i * ctxt.N;
        for (long j = 0, jdx = ctxt.Nh, idx = 0; j < slots; ++j, jdx += gap, idx += gap) {
            long mir = std::round(uvals[j].real());
            long mii = std::round(uvals[j].imag());
            mi[idx] = mir >= 0 ? (uint64_t) mir : (uint64_t) (ctxt.qVec[i] + mir);
            mi[jdx] = mii >= 0 ? (uint64_t) mii : (uint64_t) (ctxt.qVec[i] + mii);
        }
        ctxt.qiNTTAndEqual(mi, i);
    }
    delete[] uvals;
}

Plaintext my_encode_alloc(Context &ctxt, const complex<double> *v, long slots, long l) {
    auto *m = new uint64_t[l << ctxt.logN]();
    my_encode(ctxt, m, v, slots, l);
    return Plaintext(m, ctxt.N, slots, l);
}

cmpl64 *my_decrypt(Context &context, const SecretKey &secretKey, const Ciphertext &ciphertext) {
    auto msg = my_decrypt_msg(context, secretKey, ciphertext);
    return my_decode_alloc(context, msg);
}

bool isSame(const uint64_t *a, const uint64_t *b, int level, long N) {
    auto tmp = memcmp(a, b, level * N * sizeof(uint64_t));
    return tmp == 0;
}

uint64_t maxCoeffError(Context &ctxt, const uint64_t *a, const uint64_t *b, bool isNTT = true) {
    uint64_t *a_coeffs = nullptr, *b_coeffs = nullptr;
    if (isNTT) {
        a_coeffs = new uint64_t[ctxt.N];
        b_coeffs = new uint64_t[ctxt.N];
        copy(a, a + ctxt.N, a_coeffs);
        copy(b, b + ctxt.N, b_coeffs);
        ctxt.INTTAndEqual(a_coeffs, 1);
        ctxt.INTTAndEqual(b_coeffs, 1);
    } else {
        a_coeffs = const_cast<uint64_t *>(a);
        b_coeffs = const_cast<uint64_t *>(b);
    }
    uint64_t res = 0;
    for (int i = 0; i < ctxt.N; i++) {
        if (a_coeffs[i] > b_coeffs[i])
            res = std::max(res, a_coeffs[i] - b_coeffs[i]);
        else if (a_coeffs[i] < b_coeffs[i])
            res = std::max(res, b_coeffs[i] - a_coeffs[i]);
    }
    if (isNTT) {
        delete[] a_coeffs;
        delete[] b_coeffs;
    }
    return res;
}

// homomorphic evaluations
void homo_mean(Scheme &scheme, Ciphertext &ciphertext) {
    for (int i = 0; i < scheme.context.logNh; i++) {
        auto tmp = scheme.leftRotateFast(ciphertext, 1 << i);
        scheme.addAndEqual(ciphertext, tmp);
    }
    scheme.multByConstAndEqual(ciphertext, 1. / scheme.context.Nh);
    scheme.reScaleByAndEqual(ciphertext, 1);
}

Ciphertext homo_variance(Scheme &scheme, const Ciphertext &ciphertext) {
    auto mean_sqr = ciphertext, copy = ciphertext;
    scheme.conjugateAndEqual(mean_sqr); // conj(x)
    scheme.multAndEqual(mean_sqr, copy); // 0, level, ||x^2||
    scheme.reScaleByAndEqual(mean_sqr, 1); // -1 level, ||x^2||
    homo_mean(scheme, mean_sqr); // -2 level, mean(||x^2||)
    homo_mean(scheme, copy);
    auto sqr_mean = copy; // -1 level, mean(x)
    scheme.conjugateAndEqual(sqr_mean); // -1 level, conj(mean(x))
    scheme.multAndEqual(sqr_mean, copy); // -1 level, ||mean(x)^2||
    scheme.reScaleByAndEqual(sqr_mean, 1); // -2 level, ||mean(x)^2||
    scheme.subAndEqual(mean_sqr, sqr_mean); // -2 level, var(x)
    return mean_sqr;
}

Ciphertext homo_eval_poly(Scheme &scheme, const Ciphertext &ciphertext, const std::vector<double> &coeffs) {
    auto n_coeffs = coeffs.size();
    if (n_coeffs == 0)
        return ciphertext;
    if (n_coeffs == 1)
        throw std::invalid_argument("the polynomial to be evaluated should not be constant");
    auto max_deg = coeffs.size() - 1;
    auto tower_size = NTL::NumBits(max_deg);
    std::vector<Ciphertext> tower;
    tower.reserve(tower_size);
    tower.emplace_back(ciphertext);
    for (long i = 1; i < tower_size; i++) {
        Ciphertext tmp = scheme.square(tower[i - 1]);
        scheme.reScaleByAndEqual(tmp, 1);
        tower.emplace_back(tmp);
    }
    // c^(2^0), ..., c^(2^(tower_size - 1)) are computed
    Ciphertext dst = ciphertext;
    scheme.multByConstAndEqual(dst, coeffs[1]);
    scheme.reScaleByAndEqual(dst, 1);
    scheme.addConstAndEqual(dst, coeffs[0]);
    // now dst = a_0 + a_1 * x
    for (int deg = 2; deg < n_coeffs; deg++) {
        unsigned int cur_deg_total_bits = NTL::NumBits(deg), cursor_bit_idx = 0;
        for (; cursor_bit_idx < cur_deg_total_bits; cursor_bit_idx++) {
            if ((1 << cursor_bit_idx) & deg)
                break;
        }

        // we can't get the scaling factor directly, and I'm lazy to compute that
        if (fabs(coeffs[deg]) == 0) // too small s.t. encoding results is zero poly
            continue;

        Ciphertext tmp_ciphertext = tower[cursor_bit_idx];
        scheme.multByConstAndEqual(tmp_ciphertext, coeffs[deg]);
        scheme.reScaleByAndEqual(tmp_ciphertext, 1);
        while (++cursor_bit_idx < cur_deg_total_bits) {
            if ((1 << cursor_bit_idx) & deg) {
                scheme.multAndEqual(tmp_ciphertext, tower[cursor_bit_idx]);
                scheme.reScaleByAndEqual(tmp_ciphertext, 1);
            } else {
                scheme.multByConstAndEqual(tmp_ciphertext, 1);
                scheme.reScaleByAndEqual(tmp_ciphertext, 1);
            }
        }
        while (dst.l > tmp_ciphertext.l) {
            scheme.multByConstAndEqual(dst, 1);
            scheme.reScaleByAndEqual(dst, 1);
        }
        scheme.addAndEqual(dst, tmp_ciphertext);
    }
    return dst; // will RVO or NRVO optimize this out?
}

int main(int argc, char *argv[]) {
    // cmd arguments parsing
    TCLAP::CmdLine cmdLine("FullRNS-HEAAN attack&defense", ' ');
    TCLAP::ValueArg<long> seed_flag("s", "seed", "random seed", false, time(nullptr), "word-sized signed integer",
                                    cmdLine);
    TCLAP::ValueArg<uint64_t> logn_flag("n", "logn", "log_2(modulus polynomial degree)", false, 16,
                                        "positive integer <= 16", cmdLine);
    TCLAP::ValueArg<uint64_t> logp_flag("p", "logp", "log_2(scaling factor)", false, 40, "positive integer", cmdLine);
    TCLAP::ValueArg<uint64_t> levels_flag("l", "levels", "max levels", false, 10, "positive integer", cmdLine);
    TCLAP::ValuesConstraint<std::string> func_constraint({"none", "sigmoid", "exponent", "variance"});
    TCLAP::ValueArg<std::string> func_flag("f", "func", "function to evaluate", false, "none", &func_constraint,
                                           cmdLine);

    TCLAP::SwitchArg defense_flag("d", "defense", "run defense", cmdLine);

    TCLAP::SwitchArg crt_compose_flag("", "crt", "use CRT reconstruction in decryption", cmdLine);

    TCLAP::SwitchArg modified_decode_flag("", "md", "use modified decoding routine", cmdLine);

    cmdLine.parse(argc, argv);

    boolalpha(std::cout);
    auto seed = seed_flag.getValue();
    NTL::SetSeed(NTL::ZZ(seed));

    bool use_crt_construction = crt_compose_flag.getValue(), use_modified_decoding = modified_decode_flag.getValue();
    long logN = logn_flag.getValue(), logp = logp_flag.getValue(), levels = levels_flag.getValue(),
            specials = levels + 1;
    long slots = 1 << (logN - 1);
    Context context(logN, logp, levels, specials);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    scheme.addConjKey(secretKey);
    scheme.addLeftRotKeys(secretKey);

    std::cout << "total bits in ciphertext modulus: " << modulus_bits(context) << std::endl;

    std::vector<cmpl64> plain_vec;
    utils::sample_complex_circle(plain_vec, slots, 1);

    Ciphertext cipher = scheme.encrypt(plain_vec.data(), slots, levels);

    auto func = func_flag.getValue();
    utils::FuncType func_type = utils::F_NONE;
    if (func == "sigmoid")
        func_type = utils::F_SIGMOID;
    else if (func == "exponent")
        func_type = utils::F_EXPONENT;
    else if (func == "variance")
        func_type = utils::F_VARIANCE;

    // homomorphic eval
    Ciphertext homo_res;
    if (func_type == utils::F_VARIANCE)
        homo_res = homo_variance(scheme, cipher);
    else
        homo_res = homo_eval_poly(scheme, cipher, utils::func_map.at(func_type));

    std::vector<std::complex<double>> plain_vec_dst;
    double var = 0;
    if (func_type == utils::F_VARIANCE)
        utils::calc_variance(plain_vec, var);
    else
        utils::element_wise_eval_polynomial(utils::func_map.at(func_type), plain_vec,
                                            plain_vec_dst);

    // decryption
    auto dec_plain = use_crt_construction ? my_decrypt_msg(context, secretKey, homo_res) :
                     scheme.decryptMsg(secretKey, homo_res);
    auto dec_arr = use_modified_decoding ? my_decode_alloc(context, dec_plain) :
                   scheme.decode(dec_plain);

    utils::element_wise_mult_inplace(plain_vec, plain_vec);

    std::vector<cmpl64> dec_vec(dec_arr, dec_arr + slots);
    delete[] dec_arr;

    // precision estimation
    if (func_type == utils::F_VARIANCE) {
        printf("max error = %f\n", utils::max_error_1vN((std::complex<double>) var, dec_vec));
        printf("rel error = %f\n", utils::max_rel_error_1vN((std::complex<double>) var, dec_vec));
    } else {
        printf("max error = %f\n", utils::max_error(plain_vec_dst, dec_vec));
        printf("rel error = %f\n", utils::relative_error(plain_vec_dst, dec_vec));
    }


    // attack
    long common_levels = use_crt_construction ? homo_res.l : 1;
    auto re_encoded = my_encode_alloc(context, dec_vec.data(), slots, common_levels);
    std::cout << "re-encoding ok? " << isSame(re_encoded.mx, dec_plain.mx, common_levels, context.N) << '\n';

    std::cout << "max encoding error: " << maxCoeffError(context, re_encoded.mx, dec_plain.mx) << '\n';

    context.subAndEqual(re_encoded.mx, homo_res.bx, common_levels);
    ntt_inverse(context, homo_res.ax, homo_res.ax, common_levels);
    context.mulAndEqual(re_encoded.mx, homo_res.ax, common_levels);
    std::cout << "recovery ok? " << isSame(secretKey.sx, re_encoded.mx, common_levels, context.N) << '\n';
    return 0;
}