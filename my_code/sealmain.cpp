#include <seal/seal.h>
#include <vector>
#include <cstring>
#include <chrono>
#include <string>
#include "examples.h"
#include "sealattack.h"
#include "utils.h"
#include "tclap/CmdLine.h"

using namespace seal;
using namespace std;

void add_gaussian_noise(util::CoeffIter dst, const SEALContext::ContextData &ctxt, double std_dev, double max_dev){
    auto &params = ctxt.parms();
    auto coeff_count = params.poly_modulus_degree();
    auto &modulus = params.coeff_modulus();
    auto coeff_modulus_count = modulus.size();
    auto noise(util::allocate_zero_poly(coeff_count, coeff_modulus_count, MemoryManager::GetPool()));
    sealattack::sample_gaussian_noise(noise, ctxt, std_dev, max_dev);
    sealattack::add_ntt(dst, noise, dst, modulus, coeff_count);
}

double infty_norm_ntt(util::ConstCoeffIter src, const SEALContext::ContextData &ctxt){
    auto &params = ctxt.parms();
    auto coeff_count = params.poly_modulus_degree();
    auto &coeff_modulus = params.coeff_modulus();
    auto coeff_modulus_size = coeff_modulus.size();
    auto ntt_tables = ctxt.small_ntt_tables();
    auto pool = MemoryManager::GetPool();

    auto src_copy(util::allocate_zero_poly(coeff_count, coeff_modulus_size, pool));
    std::memcpy(src_copy.get(), src, coeff_count * coeff_modulus_size * sizeof(uint64_t));
    sealattack::ntt_to_crt(src_copy, coeff_modulus, coeff_count, ntt_tables);
    sealattack::crt_to_multi_precision(src_copy, ctxt, coeff_count, pool);
    std::vector<double> double_vec;
    sealattack::multi_precision_to_double_unscaled(src_copy, double_vec, ctxt);
    return utils::infty_norm(double_vec);
}

std::vector<int> modulus_bit_sizes(int total_bit_size, int special_bit_size, int scale_bit_size){
    std::vector<int> res = {special_bit_size};
    int cur_bits = special_bit_size;
    while(cur_bits < total_bit_size - special_bit_size) {
        res.push_back(scale_bit_size);
        cur_bits += scale_bit_size;
    }
    res.push_back(special_bit_size);
    return res;
}

struct OnExitTimer{
    chrono::steady_clock::time_point start;
    OnExitTimer(){
        start = chrono::steady_clock::now();
    }
    void show(const std::string& str){
        printf("%s: %f s\n", str.c_str(), chrono::duration_cast<chrono::duration<double>>(
                chrono::steady_clock::now() - start).count());
    }
    ~OnExitTimer(){
        printf("time elapsed till program ends: %f s\n", chrono::duration_cast<chrono::duration<double>>(
                chrono::steady_clock::now() - start).count());
    }
};

int main(int argc, char** argv){

    std::string mode;
    size_t noise_samp_count = 0;
    size_t logn = 0;
    size_t special_modulus_bit_length = 0;
    size_t scale_bit_length = 0;
    double std_dev = 0;
    double max_dev = 0;
    std::string infty_norm_path, mean_var_path, complete_data_path;
    try{
        TCLAP::CmdLine cmd("CKKS key recovery for Microsoft SEAL", ' ', "0.1");
        // 3 modes, attack, defense, noise
        std::vector<std::string> allowed_modes({"attack", "defense", "noise"});
        TCLAP::ValuesConstraint<std::string> mode_constraint(allowed_modes);
        TCLAP::ValueArg<std::string> mode_arg("m", "mode", "which job to run", false, "attack", &mode_constraint, cmd);

        /// for noise mode or defense mode
        // noise sample count
        TCLAP::ValueArg<std::size_t> noise_samp_count_arg("n", "samp_count", "how many noise samples are drawn", false, 1000, "a positive integer", cmd);
        // std deviation
        TCLAP::ValueArg<double> std_dev_arg("d", "std_dev", "standard deviation for gaussian distribution", false, 3.2, "a positive floating point number", cmd);
        // dist width
        TCLAP::ValueArg<std::size_t> max_dev_width_arg("w", "max_dev_width", "the ratio of max deviation and standard deviation", false, 6, "a positive integer", cmd);
        // output path for infty norms
        TCLAP::ValueArg<std::string> infty_norm_path_arg("", "inf_norm", "path to output infty norms. infty norms won't be calculated when not specified", false, "", "a path string", cmd);
        // output path for mean and variance
        TCLAP::ValueArg<std::string> mean_var_path_arg("", "mean_var", "path to output mean and variance for each slot. no calculation takes place when path is unspecified", false, "", "a path string", cmd);
        // output path for complete data
        TCLAP::ValueArg<std::string> complete_data_path_arg("", "complete_data", "path to output complete data.", false, "", "a path string", cmd);

        /// for general purposes
        // modulus polynomial degree
        TCLAP::ValueArg<std::size_t> poly_logn_arg("l", "poly_logn", "log(N) where N is the degree of cyclotomic polynomial", false, 16, "a positive integer, no greater than 16", cmd);
        // special modulus
        TCLAP::ValueArg<std::size_t> special_modulus_bit_size_arg("", "special_bits", "bit length for special modulus", false, 60, "a positive integer, no greater than 60", cmd);
        // scaling factor
        TCLAP::ValueArg<std::size_t> scale_bit_size_arg("s", "scale_bits", "bit length for scaling factor", false, 40, "a positive integer, neither greater than 60 nor special_bits", cmd);
        // parse...
        cmd.parse(argc, argv);
        // argument assignment
        mode = mode_arg.getValue();
        noise_samp_count = noise_samp_count_arg.getValue();
        logn = poly_logn_arg.getValue();
        special_modulus_bit_length = special_modulus_bit_size_arg.getValue();
        scale_bit_length = scale_bit_size_arg.getValue();
        std_dev = std_dev_arg.getValue();
        max_dev = max_dev_width_arg.getValue() * std_dev;
        infty_norm_path = infty_norm_path_arg.getValue();
        mean_var_path = mean_var_path_arg.getValue();
        complete_data_path = complete_data_path_arg.getValue();

        if(mode == "noise" && infty_norm_path.empty() && mean_var_path.empty() && complete_data_path.empty()) {
            printf("no job to do for mode \"noise\"\n");
            return 0;
        }
    }catch(TCLAP::ArgException &e){
        printf("%s for arg %s\n", e.error().c_str(), e.argId().c_str());
        return 1;
    }

    OnExitTimer timer;

    EncryptionParameters params(scheme_type::ckks);
    size_t poly_modulus_degree = 1 << logn; // value of N
    params.set_poly_modulus_degree(poly_modulus_degree);
    // ciphertext modulus q?
    params.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree,
          modulus_bit_sizes(logn == 16 ? 350 : CoeffModulus::MaxBitCount(poly_modulus_degree, sec_level_type::tc256), special_modulus_bit_length, scale_bit_length)));

    SEALContext context(params, true, sec_level_type::none);
    print_parameters(context);
    cout << endl;

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys; // evaluation key
    keygen.create_relin_keys(relin_keys);
    GaloisKeys rotation_keys;
    keygen.create_galois_keys(rotation_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    size_t slot_count = encoder.slot_count();

    vector<complex<double>> plaintext_vector;
    Plaintext plaintext_polynomial;
    double scale = size_t(1) << scale_bit_length;

    timer.show("init finished");

    if(mode == "noise"){
        auto ctxt = context.first_context_data();
//        auto &ctxt_params = ctxt->parms();
//        auto &coeff_modulus = ctxt_params.coeff_modulus();
//        auto coeff_modulus_size = coeff_modulus.size();
        auto params_id = ctxt->parms_id();

//        plaintext_polynomial.resize(poly_modulus_degree * coeff_modulus_size);
//        plaintext_polynomial.parms_id() = ctxt->parms_id();
//        plaintext_polynomial.scale() = 1;
        std::vector<double> *inf_norms = nullptr, *var = nullptr;
        std::vector<std::complex<double>> *mean = nullptr;
        FILE* complete_data_file = complete_data_path.empty() ? nullptr : fopen(complete_data_path.c_str(), "w");
        if(infty_norm_path.length())
            inf_norms = new std::vector<double>(noise_samp_count);
        if(mean_var_path.length()){
            mean = new std::vector<std::complex<double>>(slot_count);
            var = new std::vector<double>(slot_count);
        }

        std::vector<double> norm_vector;

        for(size_t i = 0; i < noise_samp_count; i++){
//            sealattack::sample_gaussian_noise(plaintext_polynomial.data(), *ctxt, 3.2, 6 * 3.2);
//            encoder.decode(plaintext_polynomial, plaintext_vector);
            encoder.fast_noise_decode(params_id, plaintext_vector, std_dev, max_dev);
            utils::to_norm(plaintext_vector, norm_vector);
            if(complete_data_file){
                for(int slot = 0; slot < slot_count; slot++)
                    fprintf(complete_data_file, "%f%s%fj ", plaintext_vector[slot].real(),
                            plaintext_vector[slot].imag() < 0 ? "" : "+", plaintext_vector[slot].imag());
                fprintf(complete_data_file, "\n");
            }
            if(inf_norms)
                (*inf_norms)[i] = utils::max_ele(norm_vector);
            if(mean){
                utils::element_wise_add_inplace(*mean, plaintext_vector); // sum up all ele
                utils::element_wise_mult_inplace(norm_vector, norm_vector);
                utils::element_wise_add_inplace(*var, norm_vector); // sum up all ||ele||^2
            }
        }

        if(complete_data_file)
            fclose(complete_data_file);
        if(inf_norms) {
            printf("mean of noise is %e\n", utils::mean(*inf_norms));
            // write to file
            auto file = fopen(infty_norm_path.c_str(), "w");
            for (auto ele : *inf_norms)
                fprintf(file, "%f\n", ele);
            fclose(file);
            delete inf_norms;
        }
        if(mean){
            utils::scale_mult(*mean, std::complex<double>(1.0/noise_samp_count));
            utils::to_norm(*mean, norm_vector);
            utils::element_wise_mult_inplace(norm_vector, norm_vector);

            utils::scale_mult(*var, 1.0/noise_samp_count);
            utils::element_wise_substract_inplace(*var, norm_vector);

            auto mean_var_file = fopen(mean_var_path.c_str(), "w");
            for(int i = 0; i < slot_count; i++){
                auto mean_ele = (*mean)[i];
                fprintf(mean_var_file, "%f%s%fj %f\n", mean_ele.real(), mean_ele.imag() < 0 ? "" : "+",
                        mean_ele.imag(), (*var)[i]);
            }
            fclose(mean_var_file);
            delete mean;
            delete var;
        }
        return 0;
    }


    printf("sampling...\n");
    utils::sample_complex_circle(plaintext_vector, slot_count, 1);
    print_vector(plaintext_vector);

    encoder.encode(plaintext_vector, scale, plaintext_polynomial);



/**
 * NOTE on encoding
 *  1.
 *   in ckks.h::encode_internal line 502, plaintext vector is encoded into a vector with DOUBLE coeff using fft
 *   in ckks.h::encode_internal line 527-619, a polynomial before NTT transformation is arranged as:
 *   (let d = coeff_modulus_size, n = coeff_count)
 *   [n uint64_t, coeffs modulo first prime] ... [...] (the count of [...] is d)
 *   >> CRT (first CRT layer) memory layout <<
 *   NTT transformation(and inverse transformation) is performed on each [n uint64_t]
 *  2.
 *  context_data.rns_tool()->base_q()->compose() transform a polynomial from CRT form into multi-precision form
 *  vice versa for decompose()
 *  3.
 *  util::ntt_negacyclic_harvey() transform a polynomial from CRT into NTT
 *  vice versa for util::inverse_ntt_negacyclic_harvey()
 *  the polynomials are by default in NTT form
 * */

    Ciphertext ciphertext;
    encryptor.encrypt(plaintext_polynomial, ciphertext);
/**
 * NOTE on encryption
 *  [Q] in encryptor.cpp line 138(encrypt_zero_internal), why is a modulus switching performed(divide_and_round_q_last_ntt_inplace)?
 *  such modulus(the last modulus in given modulus vector) is called the "special modulus", whose bit count
 *  is suggested to be set to the largest bit count of modulus chain
 *
 * NOTE on noise sampling:
 *  1. centered binomial distribution
 *   the encryption noise is sampled using centered binomial distribution (rlwe.h::sample_poly_cbd), aka CBD
 *   the CBD generates 42 random bits, subtract the sum of half of them against the other half
 *   this creates a CBD of standard variance sqrt(42/4), which is close to 3.2
 *   each coefficient is sampled according to such CBD, then transformed into CRT form,
 *   which only involves transforming negative value into upper half of modulus, since noise is far less than modulus
 *  2. discrete gaussian distribution
 *   the is not used under current configuration, which could be changed using cmake option SEAL_USE_GAUSSIAN_NOISE=ON
 *   see rlwe.h::sample_poly_normal
 *   the standard variance is the same as CBD, which is 3.2, while the max noise is 6*std_var
 *   gaussian noise is also represented in the CRT domain(first CRT layer) as the CBD noise
 * */

    evaluator.square_inplace(ciphertext);
    evaluator.relinearize_inplace(ciphertext, relin_keys);
    evaluator.rescale_to_next_inplace(ciphertext);

    Plaintext out_plaintext;
    // NOTE on the output vector: its type must be complex<double> instead of double
    //  otherwise an implicit rounding will take place, even the type of input vectors is double
    std::vector<complex<double>> out_vector;

    decryptor.decrypt(ciphertext, out_plaintext);
    encoder.decode(out_plaintext, out_vector);
    print_vector(out_vector);
    if(mode == "defense") {
        const auto &ctxt = *context.get_context_data(out_plaintext.parms_id());
        printf("scale is %e\n", out_plaintext.scale());
        printf("infty norm of plaintext is %e\n", infty_norm_ntt(out_plaintext.data(), ctxt));
        add_gaussian_noise(out_plaintext.data(), ctxt, 3.2, 6 * 3.2);
        std::vector<complex<double>> out_vector_defense;
        encoder.decode(out_plaintext, out_vector_defense);
        print_vector(out_vector_defense);
        printf("max error is %e\n", utils::max_error(out_vector, out_vector_defense));
        printf("rel error is %e\n", utils::relative_error(out_vector, out_vector_defense));
        out_vector = out_vector_defense;
    }

    // attack
    printf("start attack\n");
    auto atk_ctxt_data = context.get_context_data(out_plaintext.parms_id());
    auto coeff_modulus = atk_ctxt_data->parms().coeff_modulus();
    auto coeff_modulus_count = coeff_modulus.size();
    auto pool = MemoryManager::GetPool();
    Plaintext re_encoded;
    encoder.encode(out_vector, out_plaintext.parms_id(), out_plaintext.scale(), re_encoded);
    auto sa(util::allocate_zero_poly(poly_modulus_degree, coeff_modulus_count, pool));
    sealattack::sub_ntt(re_encoded.data(), ciphertext.data(0), sa, coeff_modulus, poly_modulus_degree);
    auto inv_a(util::allocate_zero_poly(poly_modulus_degree, coeff_modulus_count, pool));
    if(!sealattack::inv_ntt(ciphertext.data(1), inv_a, coeff_modulus, poly_modulus_degree)){
        printf("a in ciphertext in invertible\n");
        return 1;
    }
    auto recovered_sk(util::allocate_zero_poly(poly_modulus_degree, coeff_modulus_count, pool));
    sealattack::mult_ntt(inv_a, sa, recovered_sk, coeff_modulus, poly_modulus_degree);
    if(memcmp(secret_key.data().data(), recovered_sk.get(), coeff_modulus_count * poly_modulus_degree * sizeof(uint64_t)) == 0){
        printf("sk successfully recovered\n");
    }
    else{
        printf("sk recovery failed\n");
        return 1;
    }
    return 0;
}
