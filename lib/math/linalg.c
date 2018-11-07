/*
 * linalg.c
 *
 *  Created on: 05/04/2015
 *      Author: Saúl Zapotecas-Martínez
 */

#include "linalg.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "array.h"


/** ****************************************************************************
 * The following functions consider that the matrix is generated by
 * contiguous memory blocks.
 * Directives: gsl and gslcblas
 ** ***************************************************************************/

/**
 * Enforce symmetry in the matrix C
 * @param C The structure of the square matrix to enforce the symmetry
 * @param n The dimension of the square matrix
 */
void linalg_enforce_symmetry(double **C, unsigned int n)
{
    unsigned int i, j;
    double tmp;

    for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            tmp = fmax(C[i][j], C[j][i]);
            C[i][j] = C[j][i] = tmp;
        }
    }
    return;
}

/**
 * Eigen-decomposition of the matrix C
 * @param C The matrix to be decomposed. It should be positive defined.
 * @param B The normalized eigen-vectors. j-the column is the j-th eigen vector.
 * @param D The eigen-values. j-th value is the eigen-value of the j-th eigen vector.
 * @param n The dimension of C (nxn), B (nxn) and D
 */
void linalg_Eigen_decomposition(double **C, double **B, double *D, unsigned int n)
{
    unsigned int i, j;

    double *buffer;
    gsl_vector *eval = gsl_vector_alloc(n);
    gsl_matrix *evec = gsl_matrix_alloc(n, n);
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
    gsl_matrix_view C_gsl;

    buffer = (double*) malloc(sizeof(double) * (n * n));
    memcpy(buffer, C[0], sizeof(double) * (n * n));
    C_gsl = gsl_matrix_view_array(buffer, n, n);

    gsl_eigen_symmv(&C_gsl.matrix, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            B[i][j] = gsl_matrix_get(evec, i, j);
            B[j][i] = gsl_matrix_get(evec, j, i);
        }
        B[i][i] = gsl_matrix_get(evec, i, i);
        D[i] = gsl_vector_get(eval, i);
    }
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    free(buffer);
    return;
}

/**
 * Gaussian Elimination to solve Ax=b
 * Solve the linear system Ax=b (it solves for x)
 * @param A	Square Matrix (nxn)
 * @param x	The vector to be found
 * @param b	The given result
 * @param n	Dimension of x and b
 */
int linalg_Ax_b(double **A, double *x, double *b, unsigned int n)
{
    unsigned i;
    int s;
    double *buffer;
    gsl_matrix_view A_gsl;
    gsl_vector_view b_gsl;
    gsl_vector *x_gsl;
    gsl_permutation *p;
    int flag_succ;

    buffer = (double*) malloc(sizeof(double) * (n * n));
    memcpy(buffer, A[0], sizeof(double) * (n * n));

    A_gsl = gsl_matrix_view_array(buffer, n, n);
    b_gsl = gsl_vector_view_array(b, n);
    x_gsl = gsl_vector_alloc(n);
    p = gsl_permutation_alloc(n);

    gsl_linalg_LU_decomp(&A_gsl.matrix, p, &s);
    flag_succ = (fabs(gsl_linalg_LU_det(&A_gsl.matrix, s)) < 1e-10) ? 0 : 1;

    if (flag_succ == 1)
    {
        gsl_linalg_LU_solve(&A_gsl.matrix, p, &b_gsl.vector, x_gsl);
        for (i = 0; i < n; ++i)
        {
            x[i] = gsl_vector_get(x_gsl, i);
        }
    }

    /*
     printf("LU_decomp: %d\n", gsl_linalg_LU_decomp(&A_gsl.matrix, p, &s));
     printf("   LU_det: %e\n", gsl_linalg_LU_det(&A_gsl.matrix, s));
     printf(" LU_solve: %d\n",
     gsl_linalg_LU_solve(&A_gsl.matrix, p, &b_gsl.vector, x_gsl));
     */

    gsl_permutation_free(p);
    gsl_vector_free(x_gsl);

    free(buffer);
    return flag_succ;
}

/**
 * Inverse of a matrix
 * @param A 		The matrix to be inverted.
 * @param A_inv 	The inverse of C.
 * @param n     	The dimension of C (nxn), B (nxn) and D
 */
void linalg_inverse_matrix(double **A, double **A_inv, unsigned int n)
{
    unsigned int i, j;
    int s;
    double *buffer;

    gsl_matrix_view A_gsl;
    gsl_matrix *A_inv_gsl = gsl_matrix_alloc(n, n);
    gsl_permutation *perm = gsl_permutation_alloc(n);

    buffer = (double*) malloc(sizeof(double) * (n * n));
    memcpy(buffer, A[0], sizeof(double) * (n * n));
    A_gsl = gsl_matrix_view_array(buffer, n, n);

    // Make LU decomposition of matrix m
    gsl_linalg_LU_decomp(&A_gsl.matrix, perm, &s);

    // Invert the matrix m
    gsl_linalg_LU_invert(&A_gsl.matrix, perm, A_inv_gsl);

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            A_inv[i][j] = gsl_matrix_get(A_inv_gsl, i, j);
            A_inv[j][i] = gsl_matrix_get(A_inv_gsl, j, i);
        }
        A_inv[i][i] = gsl_matrix_get(A_inv_gsl, i, i);
    }

    gsl_permutation_free(perm);
    gsl_matrix_free(A_inv_gsl);
    free(buffer);
    return;
}

/**
 * Cholesky decomposition of A. A=L*L^t
 * @param A 	The square matrix to be decompose
 * @param L 	The low matrix of the decompose A
 * @param n	    The dimension of the square matrix A
 */
void linalg_Cholesky_decomposition(double **A, double **L, unsigned int n)
{
    unsigned int i, j;
    double *buffer;
    gsl_matrix_view A_gsl;

    buffer = (double*) malloc(sizeof(double) * (n * n));
    memcpy(buffer, A[0], sizeof(double) * (n * n));

    A_gsl = gsl_matrix_view_array(buffer, n, n);
    gsl_linalg_cholesky_decomp(&A_gsl.matrix);

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            L[i][j] = 0.0;
            L[j][i] = gsl_matrix_get(&A_gsl.matrix, j, i);
        }
        L[i][i] = gsl_matrix_get(&A_gsl.matrix, i, i);
    }
    free(buffer);
    return;
}

/**
 * Compute: a^{t}.cB.a
 * @param a The vector a (n)
 * @param B The square matrix B (nxn)
 * @param c The factor of B
 * @param n The dimension of vector 'a' and the square matrix 'B' (nxn)
 * @return double
 */
double linalg_at_cB_a(double *a, double **B, double c, unsigned int n)
{
    unsigned int i, j;
    double res = 0.0, sum;

    for (i = 0; i < n; ++i)
    {
        sum = 0.0; // cB x a
        for (j = 0; j < n; ++j)
        {
            sum += (c * B[i][j]) * a[j];
        }
        res += a[i] * sum; // a^{t} x cB x a
    }
    return res;
}

/**
 * Obtain the co-variance matrix 'Sigma', from the 'N' samples 'samples' of dimension 'dimension 'n'
 * @param Sigma    The resulting co-variance matrix
 * @param samples  The samples
 * @param N        The number of samples
 * @param n        The dimension of the sample
 */
void linalg_get_CM(double **Sigma, double **samples, unsigned int N, unsigned int n)
{
    unsigned int i, j, k;
    double *xk, *mu;

    /* Computing the mu vector */
    mu = (double*) malloc(sizeof(double) * n);
    for (i = 0; i < n; ++i)
    {
        mu[i] = 0.0;
        for (j = 0; j < N; ++j)
        {
            mu[i] += samples[j][i];
        }
        mu[i] = mu[i] / N;
    }

    /* Computing the Sigma matrix */
    for (i = 0; i < n; i++)
    {
        for (j = i; j < n; j++)
        {
            Sigma[i][j] = 0.0;
            for (k = 0; k < N; k++)
            {
                xk = samples[k];
                Sigma[i][j] += (xk[i] - mu[i]) * (xk[j] - mu[j]);
            }
            Sigma[i][j] = Sigma[i][j] / (N - 1.0);

            Sigma[j][i] = Sigma[i][j];
        }
    }
    free(mu);
    return;
}

/* Testing linalg_get_CM */
void Test_get_CM(void)
{
    double data[5][3] =
    {
    { 4.0, 2.0, 0.60 },
    { 4.2, 2.1, 0.59 },
    { 3.9, 2.0, 0.58 },
    { 4.3, 2.1, 0.62 },
    { 4.1, 2.2, 0.63 } };

    double **Sigma;
    double **samples;
    unsigned int n = 3U;
    unsigned int N = 5U;
    unsigned int i, j;

    Sigma = new_matrix_block(N, n);
    samples = (double**) malloc(sizeof(double*) * N);
    for (i = 0; i < N; ++i)
    {
        samples[i] = data[i];
    }

    printf("Samples:\n");
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            printf("%lf ", samples[i][j]);
        }
        printf("\n");
    }

    linalg_get_CM(Sigma, samples, N, n);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            printf("%lf ", Sigma[i][j]);
        }
        printf("\n");
    }
    free(samples);
    delete_matrix_block(Sigma);
    return;
}