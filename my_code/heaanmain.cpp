//
// Created by msh on 2021/3/10.
//

// HEAAN includes code like "using namespace std", which is really a bad habit
#include "HEAAN.h"
#include <iostream>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pE.h>
#include <NTL/GF2X.h> // ZZ_pX when p = 2
#include <NTL/GF2E.h> // ZZ_pE when p = 2
#include <ctime>
#include <exception>
#include <vector>
#include <complex>
#include <cmath>
#include "utils.h"

#define PRINTERROR(err) std::cerr << "In function " << __FUNCTION__  << " at line " << __LINE__ << ": "\
    << e.what() << std::endl;


void printParams(const Plaintext &plaintext) {
    printf("log(p) = %ld, log(q) = %ld\n", plaintext.logp, plaintext.logq);
}

/**
 * max([valid bits of coeff for coeff in arr[0:n]])
 * */
long maxBits(const NTL::ZZ *arr, long n) {
    long res = 0;
    for (long i = 0; i < n; i++) {
        res = std::max(res, NTL::NumBits(arr[i]));
    }
    return res;
}

void printPlainMag(const Plaintext &plaintext) {
    printf("magnitude of %ld\n", maxBits(plaintext.mx, plaintext.n));
}

/**
 * convert value in [0, 2^logQ - 1] to [-2^(logQ - 1), 2^(logQ - 1) - 1]
 * */
void unsignedToSigned(Plaintext &plaintext) {
    NTL::ZZ modulus = NTL::ZZ(1) << plaintext.logq;
    for (long i = 0; i < plaintext.n; i++) {
        if (NTL::NumBits(plaintext.mx[i]) == plaintext.logq)
            plaintext.mx[i] -= modulus;
    }
}

/**
 * return a complex value sampled uniformly from the unit circle
 * */
std::complex<double> sampleUnitCircle() {
    unsigned long bits = NTL::RandomWord(); // this only works on 64-bit machines
    double squareR = (double) bits / UINT64_MAX;
    bits = NTL::RandomWord();
    double theta = (double) bits / UINT64_MAX * 2 * M_PI;
    return std::polar(std::sqrt(squareR), theta);
}

/**
 * return a complex-valued array, all elements sampled uniformly from the unit circle
 * */
std::complex<double> *sampleUnitCircleArr(int n) {
    auto *vec = new std::complex<double>[n]; // bad habit, but conforms with HEAAN :)
    for (int i = 0; i < n; i++)
        vec[i] = sampleUnitCircle();
    return vec;
}

void homo_mean(Scheme &scheme, const Ciphertext &ciphertext) {
    // TODO
}

void homo_variance(Scheme &scheme, const Ciphertext &ciphertext) {
    // TODO
}

Ciphertext homo_eval_poly(Scheme &scheme, const Ciphertext &ciphertext, const std::vector<double> &coeffs) {
    auto n_coeffs = coeffs.size();
    if (n_coeffs <= 1)
        throw std::invalid_argument("the polynomial to be evaluated should not be constant");
    auto max_deg = coeffs.size() - 1;
    auto tower_size = NTL::NumBits(max_deg);
    std::vector<Ciphertext> tower;
    tower.reserve(tower_size);
    tower.emplace_back(ciphertext);
    auto log_factor = ciphertext.logp;
    for (long i = 1; i < tower_size; i++) {
        Ciphertext tmp = scheme.square(tower[i - 1]);
        scheme.reScaleByAndEqual(tmp, log_factor);
        tower.emplace_back(tmp);
    }
    // c^(2^0), ..., c^(2^(tower_size - 1)) are computed
    Ciphertext dst = ciphertext;
    scheme.multByConstAndEqual(dst, coeffs[1], log_factor);
    scheme.reScaleByAndEqual(dst, log_factor);
//    return dst;
    scheme.addConstAndEqual(dst, coeffs[0], log_factor);
    // now dst = a_0 + a_1 * x
    for (int deg = 2; deg < n_coeffs; deg++) {
        unsigned int cur_deg_total_bits = NTL::NumBits(deg), cursor_bit_idx = 0;
        for (; cursor_bit_idx < cur_deg_total_bits; cursor_bit_idx++) {
            if ((1 << cursor_bit_idx) & deg)
                break;
        }

        if (fabs(coeffs[deg]) * exp2(tower[cursor_bit_idx].logp) < 0.5) // too small s.t. encoding results is zero poly
            continue;

        Ciphertext tmp_ciphertext = tower[cursor_bit_idx];
        scheme.multByConstAndEqual(tmp_ciphertext, coeffs[deg], log_factor);
        scheme.reScaleByAndEqual(tmp_ciphertext, log_factor);
        while (++cursor_bit_idx < cur_deg_total_bits) {
            if ((1 << cursor_bit_idx) & deg) {
                scheme.multAndEqual(tmp_ciphertext, tower[cursor_bit_idx]);
                scheme.reScaleByAndEqual(tmp_ciphertext, log_factor);
            } else {
                scheme.multByConstAndEqual(tmp_ciphertext, 1, log_factor);
                scheme.reScaleByAndEqual(tmp_ciphertext, log_factor);
            }
        }
        while (dst.logq > tmp_ciphertext.logq) {
            scheme.multByConstAndEqual(dst, 1, log_factor);
            scheme.reScaleByAndEqual(dst, log_factor);
        }
        scheme.addAndEqual(dst, tmp_ciphertext);
    }
    return dst; // no move constructor, memory leak here
}

/**
 * convert ZZ-array-styled polynomials used in HEAAN into ZZX
 * */
NTL::ZZX composeZZArray(NTL::ZZ *arr, long n) {
    NTL::ZZX res;
    res.SetLength(n + 1);
    for (long i = 0; i < n; i++)
        res[i] = arr[i];
    return res; // RVO will help us
}

/**
 * convert element from ZZX(Z[x]) into ZZ_pE(R_q)
 * */
NTL::ZZ_pE convRing(const NTL::ZZX &src) {
    return NTL::conv<ZZ_pE>(NTL::conv<ZZ_pX>(src));
}

/**
 * convert element from ZZ_pE(R_q) into ZZX(Z[x])
 * */
NTL::ZZX deconvRing(const NTL::ZZ_pE &src) {
    return NTL::conv<ZZX>(NTL::conv<ZZ_pX>(src));
}

/**
 * convert element in R_q into R_2 (q is a power of 2)
 * */
NTL::GF2E convModDown(const NTL::ZZ_pE &src) {
    return NTL::conv<NTL::GF2E>(
            NTL::conv<NTL::GF2X>(
                    NTL::conv<NTL::ZZX>(
                            NTL::conv<NTL::ZZ_pX>(src)
                    )
            )
    );
}

/**
 * convert element in R_2 into R_q (q is a power of 2)
 * */
NTL::ZZ_pE convModUp(const NTL::GF2E &src) {
    return NTL::conv<NTL::ZZ_pE>(
            NTL::conv<NTL::ZZ_pX>(
                    NTL::conv<NTL::ZZX>(
                            NTL::conv<NTL::GF2X>(src)
                    )
            )
    );
}

/**
 * b + a * s = m, with a, b, s, m in R_q
 * recover s by computing a^-1 * (m - b)
 * NOTE: this requires the a is invertible in R_q
 * */
NTL::ZZ_pE recoverByInv(const NTL::ZZ_pE &m, const NTL::ZZ_pE &a, const NTL::ZZ_pE &b) {
    NTL::ZZ_pE res = m;
    res -= b;
    try {
        res /= a;
    } catch (NTL::ErrorObject &e) {
        PRINTERROR(e)
    }
    return res;
}

/**
 * b + a * s = m, with a, b, s, m in R_q
 * recover s by splitting a into a = [a] - ([a] - a), where the natural image of a in R_2
 * if [a] is invertible, [a] - a has even coefficients, i.e. [a] - a is nilpotent
 * the inverse of [a] in R_q can be constructed efficiently from its inverse in R_2
 * since [a] is invertible in R_q and [a] - a is nilpotent
 * [a] - ([a] - a) = a is also invertible in R_q
 * if we note [a] in R_q as f and ([a] - a) as g
 * then a^-1 = (f - g)^-1 = f^-1 * (1 + (g * f^-1) + (g * f^-1)^2 + ... + (g * f^-1)^(k-1)), where g^k = 0
 *
 * NOTE: this algorithm requires the (a mod 2) is invertible in R_2
 *
 * let q = 2^k, e = ceil(log_2(k))
 * this algorithm requires:
 * multiplication in R_q    2*e + 1 + 2*(e-1) + 1 + 1 = 4 * e + 1
 * addition in R_q          1 + e + 1 + (e-1) + 1 = 2 * e + 2
 * inversion in R_2         1
 *
 * the idea comes from NTRU Technical Report 9
 * https://ntru.org/f/tr/tr009v1.pdf
 * */
NTL::ZZ_pE recoverByNTRUInv(const NTL::ZZ_pE &m, const NTL::ZZ_pE &a, const NTL::ZZ_pE &b) {
    auto a2 = convModDown(a);
    NTL::GF2E a2inv;
    try {
        NTL::inv(a2inv, a2);
    } catch (NTL::ErrorObject &e) {
        PRINTERROR(e)
        return NTL::ZZ_pE();
    }
    // a = a_odd_part - a_even_part
    auto a_odd_part = convModUp(a2); // f
    auto a_even_part = a_odd_part - a; // g
    auto a_odd_part_inv = convModUp(a2inv); // f^-1 in R_2, will become f^-1 in R_q
    // compute the inverse of [a] in R_q from the inverse of [a] in R_2
    // actually computing the sequence of the inverse of a_inv_part in R_2^(2^i), where i = 0, 1, ...
    auto log2_q = NTL::NumBits(NTL::ZZ_p::modulus()) - 1;
    auto ceil_logk = NTL::NumBits(log2_q - 1); // e = ceil(log_2(k))
    for (int i = 0; i < ceil_logk; i++)
        a_odd_part_inv *= (2 - a_odd_part * a_odd_part_inv);
    // compute the inverse of a in R_q using geometric series
    auto g_finv = a_even_part * a_odd_part_inv;
    auto inv_a = g_finv + 1;
    // we don't check the exact k when g^k = 0, we use the upper bound of k = log2_q
    // let x = f^-1 * g, we then (f - g)^-1 = f^-1 * (1 + x^2 + ... + x^(k-1))
    // we also have 1 + x + ... + x^(2^e-1) = (1 + x)(1 + x^2)(1 + x^4)(1 + x^16)...(1 + x^(2^(e-1)))
    // since x^i = 0 for i >= k, we can find the smallest e satisfying 2^e - 1 >= k - 1
    // the compute the product of e = log2_k terms
    for (long cur_deg = 1; cur_deg < ceil_logk; cur_deg++) {
        // x^(2^cur_deg))
        NTL::sqr(g_finv, g_finv);
        inv_a *= g_finv + 1;
    }
    inv_a *= a_odd_part_inv;
    // inverse of a in R_q is obtained, now compute s
    return (m - b) * inv_a;
}

/**
 * b + a * s = m, with a, b, s, m in R_q
 * recover s using the same trick that Li and Micciancio have adopted
 * this trick takes advantage of the fact that sk_i = -1, 0, 1
 * NOTE: this requires the (a mod 2) is invertible in R_2
 * */
NTL::ZZ_pE recoverByTrick(const NTL::ZZ_pE &m, const NTL::ZZ_pE &a, const NTL::ZZ_pE &b) {
    auto m2 = convModDown(m), a2 = convModDown(a), b2 = convModDown(b);
    // compute [a]^-1
    NTL::GF2E a2inv;
    try {
        NTL::inv(a2inv, a2);
    } catch (NTL::ErrorObject &e) {
        PRINTERROR(e)
        return NTL::ZZ_pE();
    }
    // compute [s]
    auto s2 = (m2 - b2) * a2inv;
    // compute h, where 2 * h + [s] = s
    auto ah = convRing(deconvRing(m - b - a * convModUp(s2)) / 2);
    auto h2 = convModDown(ah) * a2inv;
    return (-2) * convModUp(h2) + convModUp(s2);
}

/**
 * perform Monte Carlo simulation to test the probability that an uniformly drawn sample from R_q is invertible
 * under the Extended Euclidean Algorithm, where q is a power of 2
 * the 'normal_bound' is the z-score of some significance level, i.e for alpha smaller but close to 1,
 * P(|x| > normal_bound) = 1 - a
 * the 'precision' is half the width of the confidence interval
 * let u = normal_bound, r = precision, the minimum sample needed 'n' will be calculated as
 *  n >= (u / 2r)^2 - u^2
 * while yielding a significance level of erf(u / sqrt(2))
 *
 * P.S. it would seem more natural to accept a parameter of significance level rather than the bound related to it
 * but the inverse of erf(or erfc) is not given in std, and I'm lazy to find one, so it comes this way
 * */
void monte_carlo_sim(double normal_bound, double precision) {
    // use the default ZZ_pE ring
    normal_bound = fabs(normal_bound);
    precision = fabs(precision);
    uint64_t n_samp = std::ceil(normal_bound * normal_bound * (1 / (4 * precision * precision) - 1));
    NTL::ZZ_pE a;
    uint64_t success_count = 0;
    for(uint64_t i = 0; i < n_samp; i++) {
        NTL::random(a);
        try {
            NTL::inv(a);
        } catch (NTL::ErrorObject& e) {
            continue;
        }
        success_count += 1;
    }
    // experiment done, compute confidence interval
    double p_observed = double(success_count) / n_samp;
//    p_observed += normal_bound * normal_bound / (2.0 * n_samp);
    double actual_precision = normal_bound * std::sqrt(p_observed * (1 - p_observed) / n_samp
            + normal_bound * normal_bound / (4.0 * n_samp * n_samp));
    double scale = n_samp / (n_samp + normal_bound * normal_bound);
    double p_guess = p_observed + normal_bound * normal_bound / (2.0 * n_samp);
    double lower = scale * (p_guess - actual_precision), upper = scale * (p_guess + actual_precision);
    double sig_level = std::erf(normal_bound / std::sqrt(2));
    printf("bound for standard normal distribution: %f, required precision: %f, significant level: %f\n"
           "sample size: %lu\n"
           "guessed EEA-invertible probability: %f\n"
           "confidence interval: (%f, %f)\n",
           normal_bound, precision, sig_level,
           n_samp,
           p_guess,
           lower, upper);
}

bool checkSame(const Plaintext &p1, const Plaintext &p2) {
    long n = p1.n;
    if (n != p2.n)
        return false;
    for (long i = 0; i < n; i++) {
        if (p1.mx[i] != p2.mx[i]) {
            std::cout << "difference at " << p1.mx[i] << " " << p2.mx[i] << std::endl;
            return false;
        }
    }
    return true;
}

int main() {
    auto timestamp = time(nullptr);
    timestamp = 1616738627;
    NTL::ZZ seed(timestamp);
    std::cout << "seed is " << seed << std::endl;
    NTL::SetSeed(seed);

    // Parameters //
    long logN = 15;
    long logQ = 353;
    long logp = 30; ///< Larger logp will give you more correct result (smaller computation noise)
    long slots = 1 << (logN - 1); ///< This should be power of two
    long numThread = 8;

    // Construct and Generate Public Keys //
    TimeUtils timeutils;
    Ring ring(logN, logQ);
    /**
     * hamming weight defaults to 64
     * each non-zero coefficient of secret key has value of 1 or -1 with uniform possibility
     * */
    SecretKey secretKey(ring);

    int total = 1, inv_success = 0, trick_success = 0, ntru_sucess = 0;

    for (int iter = 0; iter < total; iter++) {
        Scheme scheme(secretKey, ring);
        scheme.addLeftRotKeys(secretKey); ///< When you need left rotation for the vectorized message
        scheme.addRightRotKeys(secretKey); ///< When you need right rotation for the vectorized message

        // Make Random Array of Complex //
        std::complex<double> *mvec1 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);
//        std::complex<double> *mvec2 = sampleUnitCircleArr(slots);//EvaluatorUtils::randomComplexArray(slots);

        // logp seems to be the scaling factor, while logq is the ciphertext modulus
        // Encrypt Two Arry of Complex //
        Plaintext plain1 = scheme.encode(mvec1, slots, logp, logQ);
        printParams(plain1);
        printPlainMag(plain1);

        Ciphertext cipher1 = scheme.encrypt(mvec1, slots, logp, logQ);
//        Ciphertext cipher2 = scheme.encrypt(mvec2, slots, logp, logQ);

        // Addition //
//        Ciphertext cipherAdd = scheme.add(cipher1, cipher2);

        // Multiplication And Rescale //
//        Ciphertext cipherMult = scheme.mult(cipher1, cipher2);
//        Ciphertext cipherMultAfterReScale = scheme.reScaleBy(cipherMult, logp);

        // Rotation //
//        long idx = 1;
//        Ciphertext cipherRot = scheme.leftRotate(cipher1, idx);

        const std::vector<double> &coeffs = utils::func_map.at(utils::F_SIGMOID);
//        const std::vector<double> &coeffs = {0, 1};

        auto homo_res = homo_eval_poly(scheme, cipher1, coeffs);

        std::vector<std::complex<double>> plain_vec_src, plain_vec_dst;
        plain_vec_src.assign(mvec1, mvec1 + slots);
        utils::element_wise_eval_polynomial(coeffs, plain_vec_src,
                                            plain_vec_dst);

        // choose attack victim //
        Ciphertext &victim = homo_res;

        // Decrypt //
        Plaintext dmsg1 = scheme.decryptMsg(secretKey, victim);
        unsignedToSigned(dmsg1);
        std::complex<double> *dvec1 = scheme.decode(dmsg1);

        // Compare results //
        std::vector<std::complex<double>> homo_vec_res(dvec1, dvec1 + dmsg1.n);
        printf("max error = %f\n", utils::max_error(plain_vec_dst, homo_vec_res));
        printf("rel error = %f\n", utils::relative_error(plain_vec_dst, homo_vec_res));

        // attack
        // HEAAN will scale the inv-embedded polynomial by 2^(logp + ring.logQ)
        // we want to get rid of the extra scaling factor 2^(ring.logQ)
        // although round(round(m * 2^(logp + ring.logQ)) / 2^ring.logQ) will be the same as round(m * 2^logp)
        // unless the fractional part of m * 2^logp has a leading 0-bit, followed by logQ 1-bits
        // i.e m * 2^logp mod 1 = b'.01...1xxxx', where the count of bit 1 is logQ (xxx means any)
        // here is a simple trick to achieve this, without modifying library code
        Plaintext re_encoded = scheme.encode(dvec1, slots, dmsg1.logp - ring.logQ, dmsg1.logq);
        re_encoded.logp += ring.logQ;

        printParams(dmsg1);
        printPlainMag(dmsg1);
        printPlainMag(re_encoded);

        std::cout << "re-encoding ok ? " << std::boolalpha << checkSame(re_encoded, dmsg1) << std::endl;

        // init rings
        NTL::ZZ_p::init(NTL::ZZ(1) << victim.logq);
        NTL::ZZX PhiM;
        PhiM.SetLength(ring.N + 1);
        PhiM[0] = 1;
        PhiM[ring.N] = 1;
        NTL::ZZ_pE::init(NTL::conv<NTL::ZZ_pX>(PhiM));
        NTL::GF2E::init(NTL::conv<NTL::GF2X>(PhiM));

        auto m_ring = convRing(composeZZArray(re_encoded.mx, ring.N));
        auto b_ring = convRing(composeZZArray(victim.bx, ring.N));
        auto a_ring = convRing(composeZZArray(victim.ax, ring.N));

        auto recover = recoverByInv(m_ring, a_ring, b_ring);

        auto s_ring = convRing(composeZZArray(secretKey.sx, ring.N));
        std::cout << "recovery by inverse ok ? " << bool(s_ring == recover) << std::endl;

        inv_success += s_ring == recover;
        /**
         * maybe a is non-invertible, but when all of sk's coefficients are constrained into {-1, 0, 1}
         * we can found a way to decide such sk without computing a^-1
         * the method below is taken from CKKSKeyRecovery by Li and Micciancio
         *
         * 1. we know b + a * s = m, since ciphertext modulus q = 2^k,
         *  natural reduction f: R_q => R_2 forms a homomorphism (although a single-way one)
         *  ### Note: we will note the image of a polynomial 'm' in R_q under such mapping '[m]'
         *  so the following equation holds:
         *  [b] + [a] * [s] = [m]   (Eq.1)
         *  NOTE: actually the proportion of invertible elements in R_q and R_2 are the same: both are 0.5
         *   however in R_2 the inverse of any invertible element can be computed using the extended Euclidean Algorithm
         *   while the same statement doesn't hold true in R_q
         *  so with great possibility we can compute [s] = ([m] - [b]) * [a]^-1
         *  ### Note: for a polynomial 'm', we index its coefficient w.r.t x^i with 'm_i'
         *  since sk's coefficients compose of {-1, 0, 1}, [s]_i = 0 iff sk[i] = 0, and [s]_i = 1 iff sk[i] = 1 or -1
         *
         * 2. we have already distinguished all non-zero coefficients in sk, now it only remains to
         *  decide whether a non-zero coefficient is -1 or 1
         *  ### Note: we will not distinguish between [m] in R_2 and [m] in R_q, which meaning it takes depends on
         *  ###     with whom [m] is multiplied / added. e.g. in [m] * n, [m] is in R_q, while in [m] * [n], [m] is in R_2
         *  subtracting [s] from s, we have
         *  b + a * (h + [s]) = m   (Eq.2)
         *  where only h is unknown and is what we want
         *  (alternatively, we can add [s] to s, like Li and Micciancio did, which is similar)
         *  we know h_i = 0 iff s_i = 1 or 0, h_i = -2 (or 2^k - 2) iff s_i = -1
         *  all coefficients of h are even, so we can divide it by 2
         *  hx := h / 2
         *  which means
         *  b + a * (2 * hx + [s]) = m  ==>>
         *  a * hx = (m - b - a * [s]) / 2  (Eq.3)
         *  now hx_i = 0 iff s_i = 1 or 0, hx_i = -1 (or 2^(k - 1) - 1) iff s_i = -1
         *  mapping Eq.3 into R_2, we get
         *  [a] * [hx] = [(m - b - a * [s]) / 2]
         *  from which we can easily compute [hx] by multiplying both sides with [a]^-1
         *
         * 3. sk = -2 * [hx] + [s]
         *  if we defined h := s + [s], then here sk = 2 * [hx] - [s]
         *
         * 4.
         * FIXME strict proof needed
         *   a invertible polynomial 'm' in R_2 must satisfy that m(1) = 1
         *   since for any 'm' in R_2, f: R_2 -> GF(2), f(m) = m(1) is a homomorphic mapping
         *   (because 1 is a root of modulus polynomial x^n + 1)
         *   suppose m is invertible, i.e. m * m^-1 = 1, we apply f to both sides,
         *   we get m(1) * m^-1(1) = 1, then m(1) = m^-1(1) = 1
         *   what remains to be proved is that set S = {m in R_2, m(1) = 1 and m is non-invertible} is empty
         *   below gives a not-so-strict proof on this
         * Proof(sketch):
         *  consider the set S mentioned above
         *  a polynomial 'm' in R_2 is non-invertible iff gcd(x^n + 1, m) != 1
         *  while x^n + 1 = (x + 1)^n in GF(2), this can be proved by induction
         *  so if m is non-invertible, x + 1 | m
         *  we write m = (x + 1) * a, where a is some polynomial in R_2
         *  from the discussion above, we know
         *  m(1) = (x + 1)(1) * a(1) =  (1 + 1)(1) * a(1) = 0
         *  this contradicts with the assumption that m(1) = 1
         *  so the set S is empty
         *  which indicates that m is invertible iff m(1) = 1
         *  FIXME: I'm not sure if every step is correct in GF(2), for example
         *   is factorizing a polynomial in GF(2)[x] the same as that in Z[x]?
         * */
        auto recoverTrick = recoverByTrick(m_ring, a_ring, b_ring);
        std::cout << "recovery by trick ok ? " << bool(s_ring == recoverTrick) << std::endl;

        trick_success += s_ring == recoverTrick;

        auto recoverNTRU = recoverByNTRUInv(m_ring, a_ring, b_ring);
        std::cout << "recovery by NTRU-trick ok ? " << bool(s_ring == recoverNTRU) << std::endl;

        ntru_sucess += s_ring == recoverNTRU;

        delete[] mvec1;
    }

    monte_carlo_sim(3.3, 0.01);

    printf("inv success: %d\ntrick success: %d\ntotal: %d\n", inv_success, trick_success, total);

    return 0;

}