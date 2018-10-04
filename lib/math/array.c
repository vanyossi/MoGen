/*
 * array.c
 *
 *  Created on: Aug 8, 2016
 *      Author: saul
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "array.h"

/**
 * Allocate a Matrix with contiguous memory segments.
 * @param rows The number of rows of the Matrix
 * @param cols The number of columns of the Matrix
 * @return Pointer to matrix (pointer of pointers)
 */
double **new_matrix_block(unsigned int rows, unsigned int cols)
{
    unsigned int i;
    double **MAT, *buffer;

    buffer = (double*) malloc(sizeof(double) * (rows * cols));
    MAT = (double**) malloc(sizeof(double*) * rows);
    for (i = 0; i < rows; i++)
    {
        MAT[i] = &buffer[i * cols];
    }
    return MAT;
}

/**
 * Deallocate a Matrix with contiguous memory segments.
 * @param MAT Matrix pointer
 */
void delete_matrix_block(double **MAT)
{
    free(MAT[0]);
    free(MAT);
    return;
}

/**
 * Allocate memory for matrix
 * @param MAT Pointer to MAtrix
 * @param rows num of rows
 * @param cols num of columns
 */
void alloc_matrix(double **MAT, unsigned int rows, unsigned int cols)
{
    unsigned int i;

    assert(MAT == NULL);
    MAT = (double**) malloc(sizeof(double*) * rows);
    for (i = 0; i < rows; i++)
    {
        MAT[i] = (double*) malloc(sizeof(double*) * cols);
    }
    return;
}

/**
 * Free Matrix memory
 * @param MAT Matrix pointer
 * @param rows Number of rows
 */
void free_matrix(double **MAT, unsigned int rows)
{
    unsigned int i;

    for (i = 0; i < rows; i++)
    {
        free(MAT[i]);
    }
    free(MAT);
    MAT = NULL;
    return;
}

void alloc_integer_matrix(int **MAT, unsigned int rows, unsigned int cols)
{
    unsigned int i;

    assert(MAT == NULL);
    MAT = (int**) malloc(sizeof(int*) * rows);
    for (i = 0; i < rows; i++)
    {
        MAT[i] = (int*) malloc(sizeof(int*) * cols);
    }
    return;
}

void free_integer_matrix(int **MAT, unsigned int rows)
{
    unsigned int i;

    for (i = 0; i < rows; i++)
    {
        free(MAT[i]);
    }
    free(MAT);
    MAT = NULL;
    return;
}

/**
 * Set a 'value' on all the components in array 'a' (type: double)
 * @param a Pointer to array 'a'
 * @param value The value to be set
 * @param dim The dimension of the array
 */
void set_double_array(double *a, double value, int dim)
{
    int i;
    for (i = 0; i < dim; ++i)
    {
        a[i] = value;
    }
    return;
}

/**
 * Set a 'value' on all the components in array 'a' (type: integer)
 * @param a Pointer to Array 'a'
 * @param value The value to be set
 * @param dim The dimension of the array
 */
void set_integer_array(int *a, int value, int dim)
{
    int i;
    for (i = 0; i < dim; ++i)
    {
        a[i] = value;
    }
    return;
}

/**
 * Set a sequence of integer stating in 'low_bound'
 * @param a Pointer to Array 'a'
 * @param low_bound The low bound of the sequence
 * @param dim The dimension of the array
 */
void set_int_sequence(int *a, int low_bound, unsigned int dim)
{
    unsigned int i;
    int j = low_bound;
    for (i = 0; i < dim; i++)
    {
        a[i] = j++;
    }
    return;
}

/**
 * Compare array 'a' and array 'b' (type: generic)
 * @param a Array 'a'
 * @param b Array 'b'
 * @param dim The dimension of the array
 * @return 1 If 'a'=='b', 0: otherwise
 */
int equal_array(const void *a, const void *b, int dim, size_t type)
{
    int i;
    const double *d1, *d2;
    const int *i1, *i2;

    assert(dim > 0);

    if (type == sizeof(double))
    {
        d1 = (const double*) a;
        d2 = (const double*) b;
        for (i = 0; i < dim; ++i)
        {
            if (fabs(d1[i] - d2[i]) >= DBL_EPSILON)
            {
                return 0;
            }
        }
    }
    else if (type == sizeof(int))
    {
        i1 = (const int*) a;
        i2 = (const int*) b;
        for (i = 0; i < dim; ++i)
        {
            if (i1[i] != i2[i])
            {
                return 0;
            }
        }
    }
    else
    {
        printf("ERROR: No supported data (implement by yourself it)\n");
        exit(0);
    }
    return 1;
}

/**
 * Verify if a vector has a zero
 * @param fn   the vector to verified
 * @param n    the dimension of the vector
 * @return     1: if it has, 0: otherwise
 */
int has_zero(double *fn, int n)
{
    int i;
    int count = 0;
    for (i = 0; i < n; ++i)
    {
        if (fn[i] == 0.0)
        {
            count++;
            return 1;
        }
    }
    return 0;
}

/**
 * Allocate memory a vector 'V' of size 'size'
 * @param V    The vector to be allocated
 * @param size The size of the vector
 */
void alloc_INT_vector(INT_VECTOR *V, int size)
{
    V->v = (int*) malloc(sizeof(double) * (unsigned int)size);
    V->size = size;
    return;
}

/**
 * Deallocate memory a vector 'V' of size 'size'
 * @param V    The vector to be cleaned up
 */
void free_INT_vector(INT_VECTOR *V)
{
    if (V->v != NULL)
    {
        free(V->v);
    }
    V->v = NULL;
    V->size = 0;
    return;
}

void queue_integer(int *queue, int *qsize, int value)
{
    int size = (*qsize) + 1;
    queue = (int*) realloc(queue, sizeof(int) * (unsigned int)size);
    qsize[size - 1] = value;
    *qsize = size;
    return;
}

