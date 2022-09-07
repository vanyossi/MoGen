/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_GSL_VECTOR_MAP_OPS_H
#define MOGEN_GSL_VECTOR_MAP_OPS_H

#include <math.h>
#include "mgn_types.h"

double map_hpow(double value, void* pvalue);

double map_pow(double value, void* pvalue);

double map_invmul(double value, void* pvalue);

#define pmap_func_double_header(FNAME) \
double map_##FNAME(double value, void* nouse);

pmap_func_double_header(fabs)
pmap_func_double_header(sqrt)
pmap_func_double_header(acos)
pmap_func_double_header(asin)
pmap_func_double_header(atan)
pmap_func_double_header(cos)
pmap_func_double_header(sin)
pmap_func_double_header(tan)
pmap_func_double_header(acosh)
pmap_func_double_header(asinh)
pmap_func_double_header(atanh)
pmap_func_double_header(cosh)
pmap_func_double_header(sinh)
pmap_func_double_header(tanh)
pmap_func_double_header(exp)
pmap_func_double_header(exp2)
pmap_func_double_header(expm1)
pmap_func_double_header(log)
pmap_func_double_header(log10)
pmap_func_double_header(log2)
pmap_func_double_header(log1p)
pmap_func_double_header(logb)
pmap_func_double_header(ilogb)
pmap_func_double_header(cbrt)
pmap_func_double_header(erf)
pmap_func_double_header(erfc)

pmap_func_double_header(ceil)
pmap_func_double_header(floor)
pmap_func_double_header(round)
pmap_func_double_header(trunc)

//double map_fabs(double value, void* nouse) {
//    UNUSED(nouse);
//    return fabs(value);
//}

#endif //MOGEN_GSL_VECTOR_MAP_OPS_H
