//
// Created by Iván Yossi on 25/04/22.
//

#ifndef MOGEN_MGN_CEC09_H
#define MOGEN_MGN_CEC09_H

#include <gsl/gsl_vector.h>
#include <math.h>

#include "mgn_types.h"

#define cec09_unused(x) ((void)(x))


typedef struct mgn_cec09_param cec09_param;
typedef enum mgn_cec09_mop cec09_mop;

struct mgn_cec09_param {
    size_t realSize;
    size_t objSize;
    size_t consSize;
};

enum mgn_cec09_mop {
    cec_uf1 = 1 << 0,
    cec_uf2 = 1 << 1,
    cec_uf3 = 1 << 2
};

void mgn_cec09_set_limits(cec09_mop cec_mop, mgnLimit *limits)
{
    limits->min[0] = 0.0;
    limits->max[0] = 1.0;

    for (size_t i = 0; i < limits->size; ++i) {
        switch (cec_mop) {
            case cec_uf1:
            case cec_uf2:
                limits->min[i] = 0.0;
                limits->max[i] = 1.0;
                break;
            default:
                limits->min[i] = 0.0;
                limits->max[i] = 1.0;
        }
    }

}

/**
 * UF1 MOP (Uncontraint Test instances)
 * @param x Decision variable
 * @param f Objectives
 * @param g Constrains (not used)
 * @param param binary decision variable (not used)
 */
void mgn_uf1(double *x, double *f, double *g, void *param)
{
    UNUSED(g);
    cec09_param *cecp = (cec09_param*)param;

    unsigned int j, count1, count2;
    double sum1, sum2, yj;

    sum1 = sum2 = 0.0;
    count1 = count2 = 0;
    for (j = 2; j <= cecp->realSize; j++)
    {
        yj = x[j -1] - sin(6.0 * M_PI * x[0] + j * M_PI / cecp->realSize);
        yj = yj * yj;
        if (j % 2 == 0)
        {
            sum2 += yj;
            count2++;
        }
        else
        {
            sum1 += yj;
            count1++;
        }
    }
    f[0] = x[0] + 2.0 * sum1 / (double) count1;
    f[1] = 1.0 - mgn_pow(x[0],0.5) + 2.0 * sum2 / (double) count2;
}

#endif //MOGEN_MGN_CEC09_H
