/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include <stdio.h>
#include <getopt.h>

#include <math.h>
#include <gsl/gsl_matrix.h>
#include "gsl_vector_additional.h"

double cecUF123(double x)
{
    return 1 - sqrt(x);
}

double cecUF4(double x)
{
    return 1 - pow(x,2);
}

double cecUF5(double x, size_t N)
{
    return 1 - (x / (2*N));
}

double cecUF6(double x)
{
    return 1 - x;
}

double cecUF7(double x)
{
    return 1 - x;
}

int main(int argc, char const *argv[]) {
    uint16_t magic;
    uint32_t moremagic;
    size_t res = 1000;

    char ch;
    while ((ch = getopt(argc, argv, "r:")) != -1) {
        switch (ch) {
            case 'r':
                res = strtol(optarg, NULL,10);
                break;
            default:
                // nothing
                break;
        }
    }
    argc -= optind;
    argv += optind;

    double x = 0;

    gsl_matrix *m_cr = gsl_matrix_alloc(res, 2);
    res--;
    for (size_t i = 0; i <= res; ++i) {
//        gsl_matrix_set(m_cr,i,0,(double) i / (2*res));
//        gsl_matrix_set(m_cr,i,1, cecUF5(i,res));

        x = (double) i / res;
        gsl_matrix_set(m_cr,i,0,x);
        gsl_matrix_set(m_cr,i,1, cecUF7(x));
    }
     // CF6
     // intervals
//    size_t N = 2;
//    for (size_t i = 1; i <= N; ++i) {
//        double low_x = ((2.0*i) - 1) / (2 * N);
//        double high_x = (2.0 * i) / (2.0 * N);
//
//        double deltap = (high_x - low_x) / (res / N);
//        double j = low_x;
//
//        double low = ceil( (i-1)*((double)res/N));
//        double high = floor((i)*((double)res/N));
//
//        for (size_t p = low; p <= high;p++) {
////            printf("p %zu\n",p);
//            gsl_matrix_set(m_cr,p,0,j);
//            gsl_matrix_set(m_cr,p,1, cecUF6(j));
//            j += deltap;
//        }
//    }

    FILE *out = fopen("cecUF7.txt","w");
    gsl_matrix_printf(m_cr,out);
    gsl_matrix_printf(m_cr,stdout);
    fclose(out);
    gsl_matrix_free(m_cr);

    return 0;
}
