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

typedef struct point{
    size_t size;
    double *x;
} cecPoint;

void cecPoint_alloc(cecPoint *p, size_t size)
{
    p->size = size;
    p->x = calloc(size, sizeof(*p->x));
}

void cecPoint_free(cecPoint *p)
{
    free(p->x);
}

gsl_matrix* cecUF8(size_t res)
{
    res--;
    // TODO do the fast version
//        gsl_vector *xvals = gsl_vector_alloc(res);
//    for (size_t i = 0; i < res; ++i) {
//        xvals->data[i] = (double) i / res;
//    }

    // extremely inneficient FP generation
    cecPoint *points;
    printf("res: %zu %zu\n", res, sizeof(points));
    points = calloc(res,sizeof(*points));
    size_t count = 0;
    size_t count_max = res;
    printf("p: %p %zu\n", points, count_max);
    double x1, x2, x3;
    for (size_t i = 0; i <= res; ++i) {
        x1 = (double) i / res;
        for (size_t j = 0; j <= res; ++j) {
            x2 = (double) j / res;
            for (size_t k = 0; k <= res; ++k) {
                x3 = (double) k / res;
                double value = (x1 * x1) + (x2 * x2) + (x3 * x3 * x3);
                if ( fabs(value - 1) < 1e-6) {
                    cecPoint *p = &(points[count]);
                    cecPoint_alloc(p,3);
                    p->x[0] = x1;
                    p->x[1] = x2;
                    p->x[2] = x3;
                    count++;
                }
                if (count == count_max) {
                    count_max = ceil(count_max * 1.5);
                    points = realloc(points,count_max * sizeof(*points));
                }
            }
        }
    }
//    puts("");
//    printf("po: %p\n", points);
//    for (size_t i = 0; i < count; ++i) {
//        printf("%zu %p, ", i, points[i].x);
//    }
//    exit(0);
//    return NULL;
    printf("\npoints found: %zu in %p\n", count, points);
    gsl_matrix *out_FP = gsl_matrix_alloc(count,points->size);
    for (size_t i = 0; i < count; ++i) {
//        printf("%p \n", points[i].x);
        gsl_vector_view pview = gsl_vector_view_array(points[i].x, points->size);
        gsl_matrix_set_row(out_FP,i,&pview.vector);
        cecPoint_free(&points[i]);
    }
    free(points);

    return out_FP;
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


    // if cecUF08
    gsl_matrix *m_cr = cecUF8(res);
    FILE *out = fopen("cecUF8.txt","w");
    gsl_matrix_printf(m_cr,out);
    fclose(out);
    gsl_matrix_free(m_cr);
    return 0;
// uf08 stop
//
//    double x = 0;
//
//    gsl_matrix *m_cr = gsl_matrix_alloc(res, 2);
//    res--;
//    for (size_t i = 0; i <= res; ++i) {
////        gsl_matrix_set(m_cr,i,0,(double) i / (2*res));
////        gsl_matrix_set(m_cr,i,1, cecUF5(i,res));
//
//        x = (double) i / res;
//        gsl_matrix_set(m_cr,i,0,x);
//        gsl_matrix_set(m_cr,i,1, cecUF7(x));
//    }
//     // CF6
//     // intervals
////    size_t N = 2;
////    for (size_t i = 1; i <= N; ++i) {
////        double low_x = ((2.0*i) - 1) / (2 * N);
////        double high_x = (2.0 * i) / (2.0 * N);
////
////        double deltap = (high_x - low_x) / (res / N);
////        double j = low_x;
////
////        double low = ceil( (i-1)*((double)res/N));
////        double high = floor((i)*((double)res/N));
////
////        for (size_t p = low; p <= high;p++) {
//////            printf("p %zu\n",p);
////            gsl_matrix_set(m_cr,p,0,j);
////            gsl_matrix_set(m_cr,p,1, cecUF6(j));
////            j += deltap;
////        }
////    }
//
//    FILE *out = fopen("cecUF7.txt","w");
//    gsl_matrix_printf(m_cr,out);
//    gsl_matrix_printf(m_cr,stdout);
//    fclose(out);
//    gsl_matrix_free(m_cr);
//
//    return 0;
}
