//
// Created by Iván Yossi on 17/04/22.
//

#ifndef MOGEN_MGN_ZDT_H
#define MOGEN_MGN_ZDT_H

#include <gsl/gsl_vector.h>
#include <math.h>

#include "mgn_gen_operator.h" //mgn_pow

#define zdt_unused(x) ((void)(x))

//typedef enum mop_zdt_e {
//    ZDT1,
//    ZDT2,
//    ZDT3,
//    ZDT4,
//    ZDT6,
//    ZDTM1
//} ZDTVariant;

typedef struct mgn_zdt_param zdt_param;

struct mgn_zdt_param {
    size_t xsize;
    size_t fsize;
    size_t gsize;
};


void mgn_zdt1(double *x, double* f, double* g, void* sizes)
{
    zdt_unused(g);
    zdt_param *size = (zdt_param*)sizes;
    double f1, f2, gc, h, sum;

    f1 = x[0];
    sum = 0.0;
    for (size_t i = 1; i < size->xsize; i++) {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * sum / (size->xsize - 1.0);
    h = 1.0 - sqrt(f1 / gc);
    f2 = gc * h;
    f[0] = f1;
    f[1] = f2;
    return;
}

void mgn_zdt1_vector(gsl_vector *x, gsl_vector* f, gsl_vector* g, void* p)
{
    zdt_unused(g);
    zdt_unused(p);
    double f1, f2, gc, h, sum;

    f1 = fabs(gsl_vector_get(x,0));
    gsl_vector_view xall = gsl_vector_subvector(x,1,x->size-1);
    sum = gsl_vector_sum(&xall.vector);

    gc = 1.0 + 9.0 * sum / (x->size - 1.0);
    h = 1.0 - sqrt(f1 / gc);
    f2 = gc * h;
//    printf("f1, gc, h:: %.6f %.6f %.6f\n", f1, gc, h);

    gsl_vector_set(f,0,f1);
    gsl_vector_set(f,1,f2);
}

/**
 * ZDT1 Mop
 */
void mgn_zdt2(double *x, double* f, double* g, void* sizes)
{
    zdt_unused(g);
    zdt_param *size = (zdt_param*)sizes;
    double f1, f2, gc, h, sum;

    f1 = x[0];
    sum = 0.0;
    for (size_t i = 1; i < size->xsize; i++) {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * sum / (size->xsize - 1.0);
    h = 1.0 - mgn_pow(f1 / gc,2.0);
    f2 = gc * h;

    f[0] = f1;
    f[1] = f2;

    return;
}


void mgn_zdt3_vector(gsl_vector *x, gsl_vector* f, gsl_vector* g, void* p)
{
    zdt_unused(g);
    zdt_unused(p);
    double f1, f2, gc, h, sum;

    f1 = gsl_vector_get(x,0);
    gsl_vector_view xall = gsl_vector_subvector(x,1,x->size-1);
    sum = gsl_vector_sum(&xall.vector);

    gc =1.0 + 9.0 * sum / (x->size - 1.0);
    h = 1.0 - sqrt(f1 / gc) - (f1 / gc) * sin(10.0 * M_PI * f1);
    f2 = gc * h;

    gsl_vector_set(f,0,f1);
    gsl_vector_set(f,1,f2);
}

/**
 * ZDT4 Mop
 * @param mop Multiobjective problem
 * @param indv Individual to be evaluated
 */
void mgn_zdt4(double *x, double* f, double* g, void* sizes)
{
    zdt_unused(g);
    zdt_param *size = (zdt_param*)sizes;
    double f1, f2, gc, h, sum;

    f1 = x[0];
    sum = 0.0;
    for (size_t i = 1; i < size->xsize; i++)
    {
        sum += (mgn_pow(x[i], 2.0) - 10.0 * cos(4.0 * M_PI * x[i]));
    }
    gc = 1.0 + 10.0 * (size->xsize - 1.0) + sum;
    h = 1.0 - sqrt(f1 / gc);
    f2 = gc * h;

    f[0] = f1;
    f[1] = f2;
}


/**
 * ZDT6 Mop
 * @param x variable array
 * @param f solution array
 * @param g constrains array
 * @param sizes zdt_param with all array sizes
 */
void mgn_zdt6(double *x, double* f, double* g, void* sizes)
{
    zdt_unused(g);
    zdt_param *size = (zdt_param*)sizes;
    double f1, f2, gc, h, sum;

    f1 = 1.0 - exp(-4.0 * x[0]) * pow(sin(6.0 * M_PI * x[0]), 6.0);
    sum = 0.0;
    for (size_t i = 1; i < size->xsize; i++)
    {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * mgn_pow((sum / (size->xsize - 1.0)), 0.25);
    h = 1.0 - mgn_pow((f1 / gc), 2.0);
    f2 = gc * h;

    f[0] = f1;
    f[1] = f2;
}

#endif //MOGEN_MGN_ZDT_H
