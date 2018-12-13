/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *  Copyright (c) 2016 Saul Zapotecas-Martinez
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


#include "mgf_scalarization.h"

#include <stdlib.h>
#include <string.h>

#include <math.h>


void mgf_moa_set_scalarization(struct mgf_scalar_method_t *s_m, enum scalarMethod method) {
    char scalar_name[SCLM_LAST][16] = {
        "wei", "wei-norm", "weic", "weic-norm",
        "tch", "tch-norm"
    };
    scalarization_f funcs[SCLM_LAST] = {
        wei, wei_norm, weicons, weicons_norm,
        tch, tch_norm
    };

    strcpy(s_m->scalar_name, scalar_name[method]);
    s_m->func = funcs[method];
}

/**
 * @brief Weighted-Sum Approach (Gass and Satty, 1955)
 *
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double wei(int nobjs, double *f, double *w, double *z, double *n)
{
    UNUSED(z);
    UNUSED(n);
    int j;
    double sum = 0.0;
    for (j = 0; j < nobjs; j++)
    {
        sum += w[j] * f[j];
    }
    return sum;
}

/**
 * @brief Normalized Weighted-Sum Approach (Gass and Satty, 1955)
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double wei_norm(int nobjs, double *f, double *w, double *z, double *n)
{
    int j;
    double sum = 0.0;
    for (j = 0; j < nobjs; j++)
    {
        if (n[j] - z[j] == 0.0)
            return 1e+30;

        sum += w[j] * (f[j] - z[j]) / (n[j] - z[j]);
    }
    return sum;
}

/**
 * @brief Weighted-Constraint Approach (Burachik et al., 2013)
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double weicons(int nobjs, double *f, double *w, double *z, double *n)
{
    UNUSED(z);
    UNUSED(n);
    int j, i;
    double CV;
    double wifi, wjfj;
    double g;

    double sum = 0.0;

    for (i = 0; i < nobjs; i++)
    {
        wifi = w[i] * f[i];
        CV = 0.0;
        for (j = 0; j < nobjs; j++)
        {
            if (j != i)
            {
                wjfj = (w[j] * f[j]);
                g = wjfj - wifi;
                if (CV < fmax(0.0, g))
                {
                    CV = fmax(0.0, g);
                }
            }
        }
        sum += wifi + CV;
    }
    return sum;
}

/**
 * @brief Normalized Weighted Sum Contrained (Burachik et al., 2013)
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double weicons_norm(int nobjs, double *f, double *w, double *z, double *n)
{
    int j, i;
    double CV;
    double wifi, wjfj;
    double g;

    double sum = 0.0;
    double *fnorm = calloc((size_t)nobjs, sizeof(double));

    for (i = 0; i < nobjs; i++)
    {
        fnorm[i] = (f[i] - z[i]) / (n[i] - z[i]);
    }
    for (i = 0; i < nobjs; i++)
    {
        wifi = w[i] * fnorm[i];
        CV = 0.0;
        for (j = 0; j < nobjs; j++)
        {
            if (j != i)
            {
                wjfj = (w[j] * fnorm[j]);
                g = wjfj - wifi;
                if (CV < fmax(0.0, g))
                {
                    CV = fmax(0.0, g);
                }
            }
        }
        sum += wifi + CV;
    }
    free(fnorm);
    return sum;
}

/**
 * @brief Tchebycheff Approach (Bowman Jr, 1976)
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double tch(int nobjs, double *f, double *w, double *z, double *n)
{
    UNUSED(n);
    int j;
    double tmp, tch = -__DBL_MAX__;
    double wj;

    for (j = 0; j < nobjs; j++)
    {
        wj = (w[j] == 0.0) ? 0.0001 : w[j];
        tmp = wj * fabs(f[j] - z[j]);
        tch = (tch < tmp) ? tmp : tch;
    }
    return tch;
}

/**
 * @brief Normalized Tchebychef
 * @param f Objectives array
 * @param w Weights
 * @param z Max value (for nadir)
 * @param n Min value (for nadir)
 * @return weight value
 */
double tch_norm(int nobjs, double *f, double *w, double *z, double *n)
{
    int j;
    double tmp, tch = -__DBL_MAX__;
    double wj;

    for (j = 0; j < nobjs; j++)
    {
        if (n[j] - z[j] == 0.0) {
            return 1e+30;
        }
        wj = (w[j] == 0.0) ? 0.0001 : w[j];
        tmp = wj * fabs((f[j] - z[j]) / (n[j] - z[j]));
        tch = (tch < tmp) ? tmp : tch;
    }
    return tch;
}
