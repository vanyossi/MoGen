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

#include "rand.h"

#include <assert.h>
#include <string.h>

#include "math/array.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>

#include "math/gsl_extension.h"

static gsl_rng *rand_object;

/**
 * Allocate and initialize random number object/generator
 * @param seed The seed of the generator
 */
void set_random(double seed)
{
    unsigned long int uli_seed;

    assert(0.0 <= seed && seed <= 1.0);
    /* Allocate a random number generator object.
     * Different generators are available. See the docs. */
    /* second-generation version of the RANLUX algorithm of Lüscher
     * with max-decorrelation 'gsl_rng_ranlxs2', min-decorrelation 'gsl_rng_ranlxs0'
     * gsl_rng_mt19937
     * gsl_rng_taus2
     */
    if ((rand_object = gsl_rng_alloc(gsl_rng_ranlxs2)) == NULL)
    {
        printf("ERROR: Could not create random number generator\n");
        exit(1);
    }
    /* Set the seed for our generator. */
    /* The range of allowed seeds for RANLUX generator is [0,2^31-1].
     * Higher seed values are wrapped modulo 2^31.
     * Verify that: INT_MAX=2^31-1=+2147483647
     */
    uli_seed = (unsigned long int) (seed * 2147483647);
    assert(0 <= uli_seed && uli_seed <= 2147483647);
    gsl_rng_set(rand_object, uli_seed);
    return;
}

/**
 * Free up our random number
 */
void unset_random(void)
{
    gsl_rng_free(rand_object);
}

/**
 * Fetch a random double from U[0,1)
 * @return A random double from U[0,1)
 */
double rnd_perc(void)
{
    return gsl_rng_uniform(rand_object);
}

/**
 * Fetch a single random integer number in [low, high]
 * @param low lower integer bound
 * @param high high integer bound
 * @return A random integer between lower and upper
 */
int rnd_int(int low, int high)
{
    assert(low <= high);
    return (int) round(low + (rnd_perc() * (high - low)));
}

/**
 * Fetch a single random real number in [low, high)
 * @param low lower integer bound
 * @param high high integer bound
 * @return: A random real number in [low, high)
 */
double rnd_real(double low, double high)
{
    assert(low <= high);
    return (double) (low + rnd_perc() * (high - low));
}

/**
 * Fetch random number with mean and sigma
 * @param mean     The mean of the distribution
 * @param sigma    The standard deviation
 * @return A normal random number with mean and sigma
 */
double rnd_normal(double mean, double sigma)
{
    assert(sigma > 0.0);
    return (double) (mean + gsl_ran_gaussian(rand_object, sigma));
}

/**
 * Fetch a binary digit
 * @return random digit, 0 or 1
 */
int rnd_bit(void)
{
    return (rnd_perc() < 0.5) ? 1 : 0;
}

/**
 * Shuffle randomly the order of 'n' objects, each of size 'size', stored in the array 'a'
 * @param a Array of any size members
 * @param n Number of members in array
 * @param size Size of members in array.
 */
void rnd_shuffle_vector(void *a, unsigned int n, size_t size)
{
    gsl_ran_shuffle(rand_object, a, n, size);
    return;
}

/**
 * Generate a 'N' samples 'samples' from a multivaaite Gaussian distribution with covariance matrix
 * 'Sigam' and mean 'mu' of dimension 'n'
 * Multivariate Gaussian Distribution, with covariance matrix 'sigma', mean 'mu'
 * @param Sigma    Covariance matrix (nxn) (memory: necessary block gf memory)
 * @param mu       Mean
 * @param n        Dimension of the mean
 * @param samples  Samples structure (memory: don't care)
 * @param N        Number of samples
 */
void rnd_multivariate_Gaussian(double **Sigma, double *mu, unsigned int n, double **samples, unsigned int N)
{
    unsigned int i, j;
    double *buffer;
    gsl_matrix_view _Sigma;
    gsl_vector_view _mu;
    gsl_vector_view sample;

    assert(N > 0);
    assert(n > 0);

    buffer = (double*) malloc(sizeof(double) * (n * n));
    memcpy(buffer, Sigma[0], sizeof(double) * (n * n));

    _Sigma = gsl_matrix_view_array(buffer, n, n);
    _mu = gsl_vector_view_array(mu, n);

    gsl_linalg_cholesky_decomp(&_Sigma.matrix);

    /* Generating N samples */
    for (i = 0; i < N; ++i)
    {
        sample = gsl_vector_view_array(samples[i], n);
        gsl_ran_multivariate_gaussian(rand_object, &_mu.vector, &_Sigma.matrix, &sample.vector);
        for (j = 0; j < n; ++j)
        {
            samples[i][j] = gsl_vector_get(&sample.vector, j);
        }
        printf("\n");
    }
    free(buffer);
    return;
}

///**
// * Covariance calculation
// * @param r Matrix r
// * @param m Matrix m
// */
//static void cov_calculate(gsl_matrix *r, gsl_matrix *m)
//{
//    gsl_vector_view a, b;
//    size_t i, j;
//
//    for (i = 0; i < m->size1; i++)
//    {
//        for (j = 0; j < m->size2; j++)
//        {
//            double v;
//            a = gsl_matrix_column(m, i);
//            b = gsl_matrix_column(m, j);
//            v = gsl_stats_covariance(a.vector.data, a.vector.stride, b.vector.data, b.vector.stride,
//                    a.vector.size);
//            gsl_matrix_set(r, i, j, v);
//        }
//    }
//    return;
//}

/**
 * Test of the Multivariate Gaussian Distribution
 * Working well
 */
void Test_MGD()
{
    set_random(0.5);
    unsigned int n = 2;
    unsigned int N = 1000;

    double a_data[] =
    { 1, -1.5, -1.5, 3 };
    double **A;
    double *mu;
    double **samples;

    A = new_matrix_block(n, n);
    samples = new_matrix_block(N, n);
    mu = (double*) malloc(sizeof(double) * n);

    memcpy(A[0], a_data, sizeof(double) * (n * n));
    mu[0] = 0.0;
    mu[1] = 0.0;

    rnd_multivariate_Gaussian(A, mu, n, samples, N);

    unsigned int i, j;
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            printf("%lf ", samples[i][j]);
        }
        printf("\n");
    }

    delete_matrix_block(A);
    delete_matrix_block(samples);
    free(mu);
    unset_random();
    return;
}
