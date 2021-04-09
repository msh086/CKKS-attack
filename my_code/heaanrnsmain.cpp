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
    return memcmp(a, b, level * N * sizeof(uint64_t)) == 0;
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

int main() {
    boolalpha(std::cout);
    NTL::SetSeed(NTL::ZZ(12345));
    long logN = 15, logp = 40, levels = 10, specials = 11;
    long slots = 1 << (logN - 1);
    Context context(logN, logp, levels, specials);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    std::cout << "total bits in ciphertext modulus: " << modulus_bits(context) << std::endl;

    std::vector<cmpl64> plain_vec;
    utils::sample_complex_circle(plain_vec, slots, 1);

    Ciphertext cipher = scheme.encrypt(plain_vec.data(), slots, levels);
    scheme.squareAndEqual(cipher);
    scheme.reScaleByAndEqual(cipher, 1);
    auto dec_plain = my_decrypt_msg(context, secretKey, cipher);
    auto dec_arr = my_decode_alloc(context, dec_plain);
    auto dec_plain_old = scheme.decryptMsg(secretKey, cipher);
    auto dec_arr_old = scheme.decode(dec_plain_old);

    // FIXME debug
    std::cout << "CRT reconstruct unnecessary? " << isSame(dec_plain_old.mx, dec_plain.mx, 1, context.N) << '\n';
    // end FIXME

    utils::element_wise_mult_inplace(plain_vec, plain_vec);

    std::vector<cmpl64> dec_vec(dec_arr, dec_arr + slots);
    delete[] dec_arr;

    StringUtils::showcompare(plain_vec.data(), dec_vec.data(), 6, "val");

    // attack
//    auto re_encoded = scheme.encode(dec_vec.data(), slots, cipher.l);
    auto re_encoded = my_encode_alloc(context, dec_vec.data(), slots, cipher.l);
    std::cout << "re-encoding ok? " << isSame(re_encoded.mx, dec_plain.mx, cipher.l, context.N) << '\n';

    std::cout << "max encoding error: " << maxCoeffError(context, re_encoded.mx, dec_plain.mx) << '\n';

    context.subAndEqual(re_encoded.mx, cipher.bx, cipher.l);
    ntt_inverse(context, cipher.ax, cipher.ax, cipher.l);
    context.mulAndEqual(re_encoded.mx, cipher.ax, cipher.l);
    std::cout << "recovery ok? " << isSame(secretKey.sx, re_encoded.mx, cipher.l, context.N) << '\n';
    return 0;
}