/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "zdt.h"

#include <string.h>
#include <math.h>

#include "mogen_mop.h"

void mgf_zdt1(Mop *mop, Individual *indv)
{
    IndvidualType *indv_type = mgf_indv_type(indv);
    double *F = mgf_indv_get_solution_pointer(indv);
    double f1, f2, g, h, sum;

    f1 = mgf_indv_get_double(indv, 0);
    sum = 0.0;
    for (int i = 1; i < indv_type->xsize; i++) {
        sum += mgf_indv_get_double(indv, i);
    }
    g = 1.0 + 9.0 * sum / (indv_type->xsize - 1.0);
    h = 1.0 - sqrt(f1 / g);
    f2 = g * h;

    F[0] = f1;
    F[1] = f2;

    return;
}

static void mgf_zdt2(Mop *mop, Individual *indv){
    IndvidualType *indv_type = mgf_indv_type(indv);
    double *F = mgf_indv_get_solution_pointer(indv);
    double f1, f2, g, h, sum;

    f1 = mgf_indv_get_double(indv, 0);
    sum = 0.0;
    for (int i = 1; i < indv_type->xsize; i++) {
        sum += mgf_indv_get_double(indv, i);
    }
    g = 1.0 + 9.0 * sum / (indv_type->xsize - 1.0);
    h = 1.0 - pow(f1 / g, 2);
    f2 = g * h;

    F[0] = f1;
    F[1] = f2;

    return;
}

static void mgf_zdt3(Mop *mop, Individual *indv){

}

static void mgf_zdt4(Mop *mop, Individual *indv){

}

static void mgf_zdt6(Mop *mop, Individual *indv){

}

#include <stdio.h>
static void mgf_zdtm1(Mop *mop, Individual *indv)
{
    IndvidualType *indv_type = mgf_indv_type(indv);
    double *F = mgf_indv_get_solution_pointer(indv);
    double f1, f2, g, sum;

    double *real = mgf_indv_get_realdatapointer(indv);
    int *integer = mgf_indv_get_integerdatapointer(indv);

    f1 = real[0];

    sum = 0.0;
    for (int i = 1; i < indv_type->xsize; i++) {
        sum += pow(real[i] - 0.5,2);
    }
    int res;
    for (int i = 0; i < indv_type->isize; ++i) {
        res = integer[i] - (i % indv_type->isize);
        sum += res;
    }

    for (unsigned int i = 0; i < indv_type->bsize; ++i) {
        sum += mgf_indv_get_bin(indv,i) ^ (i % 2);
    }
    g = 1.0 + sum / (indv_type->xsize + indv_type->isize + indv_type->bsize - 1);
    f2 = 1 - sqrt(f1) + g;
//    g = 1.0 + 9.0 * sum / (indv_type->xsize - 1.0);
//    h = 1.0 - sqrt(f1 / g);
//    f2 = g * h;

    F[0] = f1;
    F[1] = f2;

    return;
}

struct mgf_zdt_ops_t {
    unsigned int real;
    unsigned int integer;
    unsigned int bin;
    char name[8];
    void (*zdt_eval)(Mop*, Individual*);
};

struct mgf_zdt_ops_t mgf_zdt_ops(zdt) {
    struct mgf_zdt_ops_t zdt_ops = {10,0, 0, {0}};

    if (zdt <= ZDT3) {
        zdt_ops.real = 30;
        zdt_ops.integer = 0;
        zdt_ops.bin = 0;
    }

    switch(zdt){
        case ZDT1:
            strcpy(zdt_ops.name, "ZDT1");
            zdt_ops.zdt_eval = mgf_zdt1;
            break;
        case ZDT2:
            strcpy(zdt_ops.name, "ZDT2");
            zdt_ops.zdt_eval = mgf_zdt2;
            break;
        case ZDT3:
            strcpy(zdt_ops.name, "ZDT3");
            zdt_ops.zdt_eval = mgf_zdt3;
            break;
        case ZDT4:
            strcpy(zdt_ops.name, "ZDT4");
            zdt_ops.zdt_eval = mgf_zdt4;
            break;
        case ZDT6:
            strcpy(zdt_ops.name, "ZDT6");
            zdt_ops.zdt_eval = mgf_zdt6;
            break;
        case ZDTM1:
            strcpy(zdt_ops.name, "ZDTM1");
            zdt_ops.zdt_eval = mgf_zdtm1;
            zdt_ops.real = 10;
            zdt_ops.integer = 10;
            zdt_ops.bin = 10;
            break;
    }
    return zdt_ops;
}
/**
 * The parameter settings for the unconstrained DTLZ test problems
 * @param mop Multiobjectibe prblem reference
 */
Mop *mop_zdt(ZDTVariant zdt, MopSpecs specs)
{
    struct mgf_zdt_ops_t z_ops = mgf_zdt_ops(zdt);
    Mop *mop = mogen_mop(z_ops.name, specs, 0);
    mop_set_params(mop, z_ops.real, z_ops.bin, z_ops.integer, 2, 0); // default params

    double xmin = 0.0;
    double xmax = 1.0;
    int imin = 0;
    int imax = 1;
    mop_set_limits_ndec(mop, &xmin, &xmax, 1, &imin, &imax, 1);

    mop->evaluate = z_ops.zdt_eval;

    return mop;
}
