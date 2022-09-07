/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_vector_map_ops.h"

double map_hpow(double value, void* pvalue) {
    double *p = (double*)pvalue;
    return fabs(pow(value, *p));
}

double map_pow(double value, void* pvalue) {
    double *p = (double*)pvalue;
    return pow(value, *p);
}

double map_invmul(double value, void* pvalue) {
    UNUSED(pvalue);
    return 1.0 / value;
}

#define pmap_func_double(FNAME) \
double map_##FNAME(double value, void* nouse) { \
    UNUSED(nouse); \
    return FNAME(value); \
}

pmap_func_double(fabs)
pmap_func_double(sqrt)
pmap_func_double(acos)
pmap_func_double(asin)
pmap_func_double(atan)
pmap_func_double(cos)
pmap_func_double(sin)
pmap_func_double(tan)
pmap_func_double(acosh)
pmap_func_double(asinh)
pmap_func_double(atanh)
pmap_func_double(cosh)
pmap_func_double(sinh)
pmap_func_double(tanh)
pmap_func_double(exp)
pmap_func_double(exp2)
pmap_func_double(expm1)
pmap_func_double(log)
pmap_func_double(log10)
pmap_func_double(log2)
pmap_func_double(log1p)
pmap_func_double(logb)
pmap_func_double(ilogb)
pmap_func_double(cbrt)
pmap_func_double(erf)
pmap_func_double(erfc)

pmap_func_double(ceil)
pmap_func_double(floor)
pmap_func_double(round)
pmap_func_double(trunc)
