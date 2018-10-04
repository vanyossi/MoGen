/*
 * moea_math.c
 *
 *  Created on: Aug 6, 2016
 *      Author: Saul Zapotecas-Martinez
 *     Project: Evolutionary Multi-Objective Algorithms (EMOA)
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "moea_math.h"

/**
 * The factorial of a number ('n')
 * @param n    The number to get the factorial
 * @return     The factorial of 'n'
 */
long long EMOA_factorial(long long n)
{
    // TODO, change for iterative version
    return (n <= 1LL) ? 1LL : n * EMOA_factorial(n - 1LL);
}

/**
 * The binomial coefficient (combination (k in H))
 * @param k    The number of elements to be chosen
 * @param H    The set of elements from being chosen
 * @return     The binomial coefficient (combination (k in H))
 */
unsigned int EMOA_combination(long long k, long long H)
{
    long long N, i;
    long long num = 1LL;

    for (i = 1; i <= k - 1LL; ++i)
    {
        num *= (H + i);
    }
    N = num / EMOA_factorial(k - 1LL);
    return (unsigned int) N;
}

/**
 * Euclidean Distance of two vectors ('a' and 'b')
 * @param a    The vector 'a'
 * @param b    The vector 'b0
 * @param dim  The dimension of the vectors
 * @return     The Euclidean Distance of 'a' and 'b'
 */
double EuclideanDistance(double *a, double *b, int dim)
{
    double res = 0.0;
    int i;
    for (i = 0; i < dim; i++)
    {
        res += pow(b[i] - a[i], 2.0);
    }
    return sqrt(res);
}

/**
 * p-Distance of two vectors ('a' and 'b')
 * @param a    The vector 'a'
 * @param b    The vector 'b0
 * @param dim  The dimension of the vectors
 * @param p    The p-norm to compute the distances
 * @return     The p-distance
 */
double pDistance(double *a, double *b, int dim, double p)
{
    int i;
    double dist = 0.0;

    assert(p > 0.0);
    for (i = 0; i < dim; ++i)
    {
        dist += pow(fabs(a[i] - b[i]), p);
    }
    return pow(dist, 1.0 / p);
}
