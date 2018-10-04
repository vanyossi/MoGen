/*
 *  Copyright (c) 2016 Saul Zapotecas-Martinez <saul.zapotecas@gmail.com>
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

/**
 * @brief
 * @author Saul Zapotecas-Martinez
 * @date 2106-08-05
 */

#ifndef RAND_H_
#define RAND_H_

#include <stdlib.h>

void set_random(double seed);

void unset_random(void);

double rnd_perc(void);

int rnd_int(int low, int high);

double rnd_real(double low, double high);

double rnd_normal(double mean, double sigma);

int rnd_bit(void);

void rnd_shuffle_vector(void *a, unsigned int n, size_t size);

void rnd_multivariate_Gaussian(double **Sigma, double *mu, unsigned int n, double **samples, unsigned int N);

void Test_MGD(void);

void Test_get_CM(void);

#endif /* RAND_H_ */
