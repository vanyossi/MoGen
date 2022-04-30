//
// Created by Iván Yossi on 20/04/22.
//
#include "mgn_gen_operator.h"

#include <stdio.h>

#include <math.h>
#include <float.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>

#include "mgn_random.h"

// extra
void mgn_double_exchange(double *a, double *b)
{
    double c = *a;
    *a = *b;
    *b = c;
}
// extra
// avoids nan on negative exp
double mgn_pow(double base, double nexp)
{
    double res;
    int flag = base<0;

    if (flag) {
        double rbase = (flag)? base *-1 : base;
        double rcos = cos(nexp * M_PI);

        res = pow(rbase,nexp) * ((fabs(rcos) < 1e-6 /*&& !flag*/)? 1 : rcos);
    } else {
        res = pow(base,nexp);
    }
    return res;
}

// selects 2 randon indv and gives the best one back
void mgn_select_nrandom_el(gsl_vector_int *el, gsl_vector_int *out)
{
    gsl_permutation *perm_index = gsl_permutation_alloc(out->size);
    gsl_permutation_init(perm_index);
    gsl_ran_shuffle(rnd_get_generator(),perm_index->data,perm_index->size, sizeof(size_t));
    for (size_t i = 0; i < out->size; ++i) {
        // element 0 is self
        gsl_vector_int_set(out,i, gsl_vector_int_get(
            el,gsl_permutation_get(perm_index,i))
            );
    }
    gsl_permutation_free(perm_index);
}

double fbounds(double val, double lower, double upper)
{
    return (val < lower)? lower
        : (val > upper) ? upper
        : val;
}

// Deb Agrawal 1995
//void mgn_genop_sbx(double n, double *p1, double *p2, double *c1, double *c2, size_t size
//    ,mgn_ga_sets *bounds)
//{
//    double u;
//    double Bq;
//    double nexp = 1 / (n+1);
////    printf("nexp,u,bq %.6f, %.6f, %.6f\n", nexp, u, Bq);
//
//    double cv1, cv2;
//    for (size_t i = 0; i < size; ++i) {
//        if (0.5 < rnd_getUniform()) {
//            c1[i] = p1[i];
//            c2[i] = p2[i];
//            continue;
//        }
//
//        u = rnd_getUniform();
//        Bq = (u <= 0.5)? mgn_pow(2*u,nexp) : mgn_pow(1 / (2 * (1-u)), nexp);
//
//        cv1 = 0.5 *( (1 + Bq)*p1[i] + (1 - Bq)*p2[i] );
//        cv2 = 0.5 *( (1 - Bq)*p1[i] + (1 + Bq)*p2[i] );
//
//        cv1 = fbounds(cv1,bounds->mut_llim[i], bounds->mut_ulim[i]);
//        cv2 = fbounds(cv2,bounds->mut_llim[i], bounds->mut_ulim[i]);
//
//        if (rnd_getUniform() <= 0.5) {
//            c1[i] = cv1;
//            c2[i] = cv2;
//        } else {
//            c1[i] = cv2;
//            c2[i] = cv1;
//        }
//    }
//}

void mgn_genop_sbx(double n, double *p1, double *p2, double *c1, double *c2, size_t size
                   ,mgn_ga_sets *bounds)
{
    double u;
    double Bq;
    double nexp = 1 / (n+1);
//    printf("nexp,u,bq %.6f, %.6f, %.6f\n", nexp, u, Bq);

//    double cv1, cv2;
    for (size_t i = 0; i < size; ++i) {
        if (0.5 < rnd_getUniform()) {
            c1[i] = p1[i];
            c2[i] = p2[i];
//            continue;
        } else {
            u = rnd_getUniform();
            Bq = (u <= 0.5)? mgn_pow(2*u,nexp) : mgn_pow(1 / (2 * (1-u)), nexp);

            c1[i] = 0.5 *( (1 + Bq)*p1[i] + (1 - Bq)*p2[i] );
            c2[i] = 0.5 *( (1 - Bq)*p1[i] + (1 + Bq)*p2[i] );

            c1[i] = fbounds(c1[i],bounds->mut_llim[i], bounds->mut_ulim[i]);
            c2[i] = fbounds(c2[i],bounds->mut_llim[i], bounds->mut_ulim[i]);
        }
        if (rnd_getUniform() > 0.5) {
            mgn_double_exchange(&c1[i], &c2[i]);
        }
    }
}

void mgn_genop_sbx_alt(double eta_c, double *par1, double *par2, double *chi1, double *chi2, size_t size
                       ,mgn_ga_sets *params)
{
    double *low_bound = params->mut_llim;
    double *up_bound = params->mut_ulim;

    unsigned int i;
    double rand;
    double y1, y2, yl, yu;
    double cv1, cv2;
    double alpha, beta, betaq;

//    if (rnd_getUniform() > Pc)
//    {
//        for (i = 0; i < size; i++)
//        {
//            chi1[i] = par1[i];
//            chi2[i] = par2[i];
//        }
//        return;
//    }
    for (i = 0; i < size; i++) {
//        if (0.5 < rnd_getUniform()) {
//            chi1[i] = par1[i];
//            chi2[i] = par2[i];
//            continue;
//        };

        if (rnd_getUniform() <= 0.5) {
            if (fabs(par1[i] - par2[i]) > 1.0e-14) {
                if (par1[i] < par2[i]) {
                    y1 = par1[i];
                    y2 = par2[i];
                } else {
                    y1 = par2[i];
                    y2 = par1[i];
                }
                yl = low_bound[i];
                yu = up_bound[i];
                rand = rnd_getUniform() + .000001;
                beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                alpha = 2.0 - mgn_pow(beta, -(eta_c + 1.0));
                if (rand <= (1.0 / alpha)) {
                    betaq = mgn_pow((rand * alpha), (1.0 / (eta_c + 1.0)));
                } else {
                    betaq = mgn_pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                }
                cv1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                alpha = 2.0 - mgn_pow(beta, -(eta_c + 1.0));
                if (rand <= (1.0 / alpha)) {
                    betaq = mgn_pow((rand * alpha), (1.0 / (eta_c + 1.0)));
                } else {
                    betaq = mgn_pow((1.0 / (2.0 - rand * alpha)), (1.0 / (eta_c + 1.0)));
                }
                cv2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

                if (cv1 < yl)
                    cv1 = yl;
                if (cv2 < yl)
                    cv2 = yl;
                if (cv1 > yu)
                    cv1 = yu;
                if (cv2 > yu)
                    cv2 = yu;
                if (rnd_getUniform() <= 0.5) {
                    chi1[i] = cv2;
                    chi2[i] = cv1;
                } else {
                    chi1[i] = cv1;
                    chi2[i] = cv2;
                }
            } else {
                chi1[i] = par1[i];
                chi2[i] = par2[i];
            }
        } else {
            chi1[i] = par1[i];
            chi2[i] = par2[i];
        }
    }
}

/*
 *
 *    Mutation
 *
 */

// deb 2014
void mgn_genop_pbm(double n, double pm, double *p, const double *lb, const double *ub, size_t size)
{
    double nexp = 1 / (n+1);
    double u;
    double np;
    double d;

    // TODO mutate for each
    for (size_t i = 0; i < size; ++i) {
        if (pm < rnd_getUniform()) continue;

        u = rnd_getUniform();
        if (u <= 0.5) {
            d = mgn_pow((2*u),nexp)-1;
            d = (isnan(d))? DBL_MAX : d;
            np = p[i] + d * (p[i] - lb[i]);
        } else {
            d = 1-(2 * (1-u)) * nexp;
            np = p[i] + d * (ub[i] - p[i]);
        }
        np = (np < lb[i])? lb[i]: np;
        np = (np > ub[i])? ub[i]: np;
        p[i] = np;
    }
}

/*
void PBM(INDIVIDUAL *ind, double *low_bound, double *up_bound, double pbm_eta, double Pm)
{
    unsigned int j;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (j = 0; j < mop.nreal; j++)
    {
        if (rnd_perc() <= Pm)
        {
            y = ind->x[j];
            yl = low_bound[j];
            yu = up_bound[j];
            delta1 = (y - yl) / (yu - yl);
            delta2 = (yu - y) / (yu - yl);
            rnd = rnd_perc();
            mut_pow = 1.0 / (pbm_eta + 1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0 - delta1;
                val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (pbm_eta + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0 - delta2;
                val = 2.0 * (1.0 - rnd)
                      + 2.0 * (rnd - 0.5) * (pow(xy, (pbm_eta + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }
            y = y + deltaq * (yu - yl);
            if (y < yl)
                y = yl;
            if (y > yu)
                y = yu;
            ind->x[j] = y;
        }
    }
    return;
}
*/

