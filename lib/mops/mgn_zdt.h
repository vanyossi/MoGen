//
// Created by Iván Yossi on 17/04/22.
//

#ifndef MOGEN_MGN_ZDT_H
#define MOGEN_MGN_ZDT_H

#include <math.h>
#include <string.h>
#include <gsl/gsl_vector.h>

#include "mgn_mop.h"
#include "mgn_gen_operator.h" //mgn_pow

#define zdt_unused(x) ((void)(x))

typedef struct mgn_zdt_param zdt_param;

struct mgn_zdt_param {
    mgn_mop_param_common()
    size_t x_size;
    size_t f_size;
    size_t g_size;
};


void mgn_zdt1(double *x, double* f, double* g, void* sizes)
{
    zdt_unused(g);
    zdt_param *size = (zdt_param*)sizes;
    double f1, f2, gc, h, sum;

    f1 = x[0];
    sum = 0.0;
    for (size_t i = 1; i < size->x_size; i++) {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * sum / (size->x_size - 1.0);
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

//    for (size_t i = 0; i < x->size; ++i) {
//        printf("%g, ", gsl_vector_get(x,i));
//    }
//    puts("");
//    printf("f1 f2, gc, h:: %.6f %.6f %.6f\n", f1,f2, gc, h);

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
    for (size_t i = 1; i < size->x_size; i++) {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * sum / (size->x_size - 1.0);
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
    for (size_t i = 1; i < size->x_size; i++)
    {
        sum += (mgn_pow(x[i], 2.0) - 10.0 * cos(4.0 * M_PI * x[i]));
    }
    gc = 1.0 + 10.0 * (size->x_size - 1.0) + sum;
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
    for (size_t i = 1; i < size->x_size; i++)
    {
        sum += x[i];
    }
    gc = 1.0 + 9.0 * mgn_pow((sum / (size->x_size - 1.0)), 0.25);
    h = 1.0 - mgn_pow((f1 / gc), 2.0);
    f2 = gc * h;

    f[0] = f1;
    f[1] = f2;
}

typedef enum mop_zdt_e {
    ZDT1,
    ZDT2,
    ZDT3,
    ZDT4,
    ZDT5,
    ZDT6,
    ZDTM1
} MGN_ZDT_VAR;


MGN_ZDT_VAR mop_zdt_str_toenum(char *type)
{

    MGN_ZDT_VAR etype = 0;
    if ( strcasecmp(type,"ZDT1") == 0 ) {
        etype = ZDT1;
    } else if ( strcasecmp(type,"ZDT2") == 0 ) {
        etype = ZDT2;
    } else if ( strcasecmp(type,"ZDT3") == 0 ) {
        etype = ZDT3;
    } else if ( strcasecmp(type,"ZDT4") == 0 ) {
        etype = ZDT4;
    } else if ( strcasecmp(type,"ZDT5") == 0 ) {
        etype = ZDT5;
    } else if ( strcasecmp(type,"ZDT6") == 0 ) {
        etype = ZDT6;
    } else if ( strcasecmp(type,"ZDTM1") == 0 ) {
        etype = ZDTM1;
    } else {
        etype = ZDT1;
    }
    return etype;
}

void mgn_zdt_free(mgnMop* mop)
{
    zdt_param* cp = (zdt_param*)mop->params;
    mgn_limit_free(mop->limits);
    free(cp);
}

mgnLimit* mgn_zdt_alloc_limits(MGN_ZDT_VAR variant, size_t x_size)
{
    UNUSED(variant);
    mgnLimit *moplim = mgn_limit_alloc(x_size);
    for (size_t i = 0; i < moplim->size; ++i) {
        moplim->min[i] = 0;
        moplim->max[i] = 1;
    }
    return moplim;
}

mgnMop* mgn_zdt_init(MGN_ZDT_VAR variant, mgn_indv_param *param)
{
    mgnMop *mop = mgn_mop_alloc(param);
    mop->free = mgn_zdt_free;
    zdt_param* cp = malloc(sizeof(*cp));
    mop->limits = mgn_zdt_alloc_limits(variant, param->x_size);
    cp->pos = 0;
    cp->x_size = param->x_size;
    cp->f_size = param->f_size;
    cp->g_size = param->g_size;

    mop->params = cp;

    switch (variant) {
        case ZDT1:
            strcpy(mop->name, "ZDT1");
            mop->eval_array = mgn_cast_eval(mgn_zdt1);
            mop->eval = mgn_cast_eval(mgn_zdt1_vector);
            break;
        case ZDT2:
            strcpy(mop->name, "ZDT2");
            mop->eval_array = mgn_cast_eval(mgn_zdt2);
//            mop->eval = mgn_cast_eval(mgn_zdt2_vector);
            break;
        case ZDT3:
            strcpy(mop->name, "ZDT3");
//            mop->eval_array = mgn_cast_eval(mgn_zdt3);
            mop->eval = mgn_cast_eval(mgn_zdt3_vector);
            break;
        case ZDT4:
            strcpy(mop->name, "ZDT4");
            mop->eval_array = mgn_cast_eval(mgn_zdt4);
//            mop->eval = mgn_cast_eval(mgn_zdt1_vector);
            break;
        case ZDT6:
            strcpy(mop->name, "ZDT6");
            mop->eval_array = mgn_cast_eval(mgn_zdt6);
//            mop->eval = mgn_cast_eval(mgn_zdt1_vector);
            break;
        default:
            break;
    }
    return mop;
}

#endif //MOGEN_MGN_ZDT_H
