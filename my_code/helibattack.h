//
// Created by msh on 2021/3/2.
//

#ifndef CKKS_HELIBATTACK_H
#define CKKS_HELIBATTACK_H

#include <helib/helib.h>
#include <NTL/xdouble.h>

namespace helibattack {
    void eval_polynomial(const std::vector<double>& coeffs, const helib::Ctxt& x, helib::Ctxt& dst);

    void eval_polynomial_inplace(const std::vector<double>& coeffs, helib::Ctxt& x);

    // FIXME: for debug only, remove me later
    void print_NTL_macros();

    void print_xdouble(NTL::xdouble val);
}


#endif //CKKS_HELIBATTACK_H
