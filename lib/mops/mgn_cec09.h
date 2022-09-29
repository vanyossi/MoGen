//
// Created by Iván Yossi on 25/04/22.
//

#ifndef MOGEN_MGN_CEC09_H
#define MOGEN_MGN_CEC09_H

#include <gsl/gsl_vector.h>
#include <math.h>
#include <string.h>

#include "mgn_types.h"
#include "mgn_gen_operator.h"

//#define cec09_unused(x) ((void)(x))

typedef struct mgn_cec09_param cec09_param;


// TODO finish migration to param common
struct mgn_cec09_param {
    mgn_mop_param_common()
    size_t x_size;
    size_t f_size;
    size_t g_size;
};

typedef enum mop_cec09_e {
    UF1 = 1 << 0,
    UF2,
    UF3,
    UF4,
    UF5,
    UF6,
    UF7,
    UF8,
    UF9,
    UF10
} MGN_CEC09_VAR;

MGN_CEC09_VAR mop_cec09_str_toenum(char *type)
{

    MGN_CEC09_VAR etype = 0;
    if ( strcasecmp(type,"UF1") == 0 ) {
        etype = UF1;
    } else if ( strcasecmp(type,"UF2") == 0 ) {
        etype = UF2;
    } else if ( strcasecmp(type,"UF3") == 0 ) {
        etype = UF3;
    } else if ( strcasecmp(type,"UF4") == 0 ) {
        etype = UF4;
    } else if ( strcasecmp(type,"UF5") == 0 ) {
        etype = UF5;
    } else if ( strcasecmp(type,"UF6") == 0 ) {
        etype = UF6;
    } else if ( strcasecmp(type,"UF7") == 0 ) {
        etype = UF7;
    } else if ( strcasecmp(type,"UF8") == 0 ) {
        etype = UF8;
    } else if ( strcasecmp(type,"UF9") == 0 ) {
        etype = UF9;
    } else {
        etype = UF10;
    }
    return etype;
}

void mgn_cec09_set_limits(MGN_CEC09_VAR cec_mop, mgnLimit *limits)
{
    limits->min[0] = 0.0;
    limits->max[0] = 1.0;

    for (size_t i = 1; i < limits->size; ++i) {
        switch (cec_mop) {
            case UF1:
                limits->min[i] = 0.0;
                limits->max[i] = 1.0;
            case UF2:
                limits->min[i] = -1.0;
                limits->max[i] = 1.0;
                break;
            case UF3:
            case UF4:
            case UF5:
            case UF6:
            case UF7:
            case UF9:
            case UF10:
                limits->min[i] = -2.0;
                limits->max[i] = 2.0;
                break;
            case UF8:
                limits->min[i] = -4.0;
                limits->max[i] = 4.0;
                break;
            default:
                limits->min[i] = 0.0;
                limits->max[i] = 0.0;
        }
    }

    if (cec_mop == UF8 || cec_mop == UF9 || cec_mop == UF10) {
        limits->min[1] = 0.0;
        limits->max[1] = 1.0;
    }
    // TODO missing contrains
}


/**
 * UF1 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param param binary decision variable (not used)
 */
void mgn_uf1(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;

    unsigned int j, count1, count2;
    double sum1, sum2, yj;

    sum1 = sum2 = 0.0;
    count1 = 0;
    count2 = 0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j -1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        yj = yj * yj;

        if (j % 2 == 0){
            sum2 += yj;
            count2++;
        }
        else {
            sum1 += yj;
            count1++;
        }
    }

    f[0] = x[0] + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - mgn_pow(x[0],0.5) + 2.0 * sum2 / (double) count2;
}

/**
 * UF2 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf2(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;

    double sum1, sum2, yj;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        if (j % 2 == 0)
        {
            yj = x[j - 1]
                 - 0.3 * x[0]
                   * (x[0] * cos(24.0 * M_PI * x[0] + 4.0 * j * M_PI / cecp->x_size) + 2.0)
                   * sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            yj = x[j - 1]
                 - 0.3 * x[0]
                   * (x[0] * cos(24.0 * M_PI * x[0] + 4.0 * j * M_PI / cecp->x_size) + 2.0)
                   * cos(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
            sum1 += yj * yj;
            count1++;
        }
    }
    f[0] = x[0] + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double) count2;
    return;
}

/**
 * UF3 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf3(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;
    double sum1, sum2, prod1, prod2, yj, pj;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    prod1 = prod2 = 1.0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - pow(x[0], 0.5 * (1.0 + 3.0 * (j - 2.0) / (cecp->x_size - 2.0)));
        pj = cos(20.0 * yj * M_PI / sqrt(j + 0.0));
        if (j % 2 == 0)
        {
            sum2 += yj * yj;
            prod2 *= pj;
            count2++;
        }
        else
        {
            sum1 += yj * yj;
            prod1 *= pj;
            count1++;
        }
    }
    f[0] = x[0] + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double) count1;
    f[1] = 1.0 - sqrt(x[0]) + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double) count2;
    return;
}

/**
 * UF4 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf4(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;
    double sum1, sum2, yj, hj;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        hj = fabs(yj) / (1.0 + exp(2.0 * fabs(yj)));
        if (j % 2 == 0)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum1 += hj;
            count1++;
        }
    }
    f[0] = x[0] + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - x[0] * x[0] + 2.0 * sum2 / (double) count2;
    return;
}

/**
 * UF5 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf5(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;
    double sum1, sum2, yj, hj, N, E;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    N = 10.0;
    E = 0.1;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        hj = 2.0 * yj * yj - cos(4.0 * M_PI * yj) + 1.0;
        if (j % 2 == 0)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum1 += hj;
            count1++;
        }
    }
    hj = (0.5 / N + E) * fabs(sin(2.0 * N * M_PI * x[0]));
    f[0] = x[0] + hj + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - x[0] + hj + 2.0 * sum2 / (double) count2;
    return;
}

/**
 * UF6 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf6(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;
    double sum1, sum2, prod1, prod2, yj, hj, pj, N, E;

    N = 2.0;
    E = 0.1;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    prod1 = prod2 = 1.0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        pj = cos(20.0 * yj * M_PI / sqrt(j + 0.0));
        if (j % 2 == 0)
        {
            sum2 += yj * yj;
            prod2 *= pj;
            count2++;
        }
        else
        {
            sum1 += yj * yj;
            prod1 *= pj;
            count1++;
        }
    }

    hj = 2.0 * (0.5 / N + E) * sin(2.0 * N * M_PI * x[0]);
    if (hj < 0.0)
        hj = 0.0;
    f[0] = x[0] + hj + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double) count1;
    f[1] = 1.0 - x[0] + hj + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double) count2;
    return;
}

/**
 * UF7 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf7(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2;
    double sum1, sum2, yj;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    for (j = 2; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        if (j % 2 == 0)
        {
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            sum1 += yj * yj;
            count1++;
        }
    }
    yj = pow(x[0], 0.2);
    f[0] = yj + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - yj + 2.0 * sum2 / (double) count2;
    return;
}

/**
 * UF8 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf8(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2, count3;
    double sum1, sum2, sum3, yj;

    sum1 = sum2 = sum3 = 0.0;
    count1 = count2 = count3 = 0;
    for (j = 3; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - 2.0 * x[1] * sin(2.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        if (j % 3 == 1)
        {
            sum1 += yj * yj;
            count1++;
        }
        else if (j % 3 == 2)
        {
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            sum3 += yj * yj;
            count3++;
        }
    }
    f[0] = cos(0.5 * M_PI * x[0]) * cos(0.5 * M_PI * x[1]) + 2.0 * sum1 / (double) count1;
    f[1] = cos(0.5 * M_PI * x[0]) * sin(0.5 * M_PI * x[1]) + 2.0 * sum2 / (double) count2;
    f[2] = sin(0.5 * M_PI * x[0]) + 2.0 * sum3 / (double) count3;
    return;
}

/**
 * UF9 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf9(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2, count3;
    double sum1, sum2, sum3, yj, E;

    E = 0.1;
    sum1 = sum2 = sum3 = 0.0;
    count1 = count2 = count3 = 0;
    for (j = 3; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - 2.0 * x[1] * sin(2.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        if (j % 3 == 1)
        {
            sum1 += yj * yj;
            count1++;
        }
        else if (j % 3 == 2)
        {
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            sum3 += yj * yj;
            count3++;
        }
    }
    yj = (1.0 + E) * (1.0 - 4.0 * (2.0 * x[0] - 1.0) * (2.0 * x[0] - 1.0));
    if (yj < 0.0)
        yj = 0.0;
    f[0] = 0.5 * (yj + 2 * x[0]) * x[1] + 2.0 * sum1 / (double) count1;
    f[1] = 0.5 * (yj - 2 * x[0] + 2.0) * x[1] + 2.0 * sum2 / (double) count2;
    f[2] = 1.0 - x[1] + 2.0 * sum3 / (double) count3;
    return;
}

/**
 * UF10 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param params binary decision variable (not used)
 */
void mgn_uf10(double *x, double *f, double *g, void *params)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)params;
    unsigned int j, count1, count2, count3;
    double sum1, sum2, sum3, yj, hj;

    sum1 = sum2 = sum3 = 0.0;
    count1 = count2 = count3 = 0;
    for (j = 3; j <= cecp->x_size; j++)
    {
        yj = x[j - 1] - 2.0 * x[1] * sin(2.0 * M_PI * x[0] + j * M_PI / cecp->x_size);
        hj = 4.0 * yj * yj - cos(8.0 * M_PI * yj) + 1.0;
        if (j % 3 == 1)
        {
            sum1 += hj;
            count1++;
        }
        else if (j % 3 == 2)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum3 += hj;
            count3++;
        }
    }
    f[0] = cos(0.5 * M_PI * x[0]) * cos(0.5 * M_PI * x[1]) + 2.0 * sum1 / (double) count1;
    f[1] = cos(0.5 * M_PI * x[0]) * sin(0.5 * M_PI * x[1]) + 2.0 * sum2 / (double) count2;
    f[2] = sin(0.5 * M_PI * x[0]) + 2.0 * sum3 / (double) count3;
    return;
}


void pmgn_cec09_free(mgnMop* mop)
{
    cec09_param* cp = (cec09_param*)mop->params;
    mgn_limit_free(mop->limits);
    free(cp);
}

// TODO boundary settings (?)
mgnMop* mgn_cec09_init(MGN_CEC09_VAR variant, mgn_indv_param * param)
{
    mgnMop *mop = mgn_mop_alloc(param);
    mop->free = pmgn_cec09_free;

    mgnLimit *limit = mgn_limit_alloc(param->x_size);
    mgn_cec09_set_limits(variant, limit);
    mop->limits = limit;

    cec09_param* cp = malloc(sizeof(*cp));

    cp->pos = 0;
    cp->x_size = param->x_size;
    cp->f_size = param->f_size;
    cp->g_size = param->g_size;

    mop->params = cp;

    switch (variant) {
        case UF1:
            strcpy(mop->name, "UF1");
            mop->eval_array = mgn_cast_eval(mgn_uf1);
//            mop->eval = mgn_cast_eval(mgn_uf1_vector);
            break;

        case UF2:
            strcpy(mop->name, "UF2");
            mop->eval_array = mgn_cast_eval(mgn_uf2);
            break;

        case UF3:
            strcpy(mop->name, "UF3");
            mop->eval_array = mgn_cast_eval(mgn_uf3);
            break;

        case UF4:
            strcpy(mop->name, "UF4");
            mop->eval_array = mgn_cast_eval(mgn_uf4);
            break;

        case UF5:
            strcpy(mop->name, "UF5");
            mop->eval_array = mgn_cast_eval(mgn_uf5);
            break;

        case UF6:
            strcpy(mop->name, "UF6");
            mop->eval_array = mgn_cast_eval(mgn_uf6);
            break;

        case UF7:
            strcpy(mop->name, "UF7");
            mop->eval_array = mgn_cast_eval(mgn_uf7);
            break;

        case UF8:
            strcpy(mop->name, "UF8");
            mop->eval_array = mgn_cast_eval(mgn_uf8);
            break;

        case UF9:
            strcpy(mop->name, "UF9");
            mop->eval_array = mgn_cast_eval(mgn_uf9);
            break;

        case UF10:
            strcpy(mop->name, "UF10");
            mop->eval_array = mgn_cast_eval(mgn_uf10);
            break;

        default:
            break;
    }
    return mop;
}

#endif //MOGEN_MGN_CEC09_H
