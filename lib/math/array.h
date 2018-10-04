/*
 * array.h
 *
 *  Created on: Aug 8, 2016
 *      Author: saul
 */

#ifndef MATH_ARRAY_H_
#define MATH_ARRAY_H_

typedef struct
{
    int *v;
    int size;
} INT_VECTOR;

double **new_matrix_block(unsigned int rows, unsigned int cols);
void delete_matrix_block(double **MAT);

void alloc_matrix(double **MAT, unsigned int rows, unsigned int cols);
void free_matrix(double **MAT, unsigned int rows);
void alloc_integer_matrix(int **MAT, unsigned int rows, unsigned int cols);
void free_integer_matrix(int **MAT, unsigned int rows);


void set_double_array(double *a, double value, int dim);
void set_integer_array(int *a, int value, int dim);
void set_int_sequence(int *a, int low_bound, unsigned int dim);
int equal_array(const void *a, const void *b, int dim, size_t size);

int has_zero(double *fn, int n);
void alloc_INT_vector(INT_VECTOR *V, int size);
void free_INT_vector(INT_VECTOR *V);
void queue_integer(int *queue, int *qsize, int value);


#endif /* MATH_ARRAY_H_ */
