/*
 * math.h
 *
 *  Created on: Aug 6, 2016
 *      Author: Saul Zapotecas-Martinez
 *     Project: Evolutionary Multi-Objective Algorithms (EMOA)
 */

#ifndef MATH_MATH_H_
#define MATH_MATH_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long long EMOA_factorial(long long n);
unsigned int EMOA_combination(long long k, long long H);
double EuclideanDistance(double *a, double *b, int dim);
double pDistance(double *a, double *b, int dim, double p);


#endif /* MATH_MATH_H_ */
