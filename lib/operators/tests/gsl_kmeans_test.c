/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_kmeans_new.h"
#include "mgn_random.h"

#include <string.h>

void gsl_matrix_printf(gsl_matrix *M, FILE *stream)
{
    const char *format = "%.6f";

    for (int c = 0; c < M->size1; ++c) {
        gsl_vector_view crow = gsl_matrix_row(M, c);
        for (size_t val = 0; val < crow.vector.size; ++val) {
            fprintf(stream,format, gsl_vector_get(&crow.vector, val));
            if(val != crow.vector.size -1) {
                fprintf(stream," ");
            }
        }
        fprintf(stream, "\n");
    }
}

int main(int argc, char const *argv[]) {

    rnd_initialize();
    rnd_set_seed(23);

    gsl_matrix *X = gsl_matrix_alloc(6,2);

    gsl_matrix_set(X,0,0,1);
    gsl_matrix_set(X,0,1,6);

    gsl_matrix_set(X,1,0,4);
    gsl_matrix_set(X,1,1,8);

    gsl_matrix_set(X,2,0,1);
    gsl_matrix_set(X,2,1,9);

    gsl_matrix_set(X,3,0,6);
    gsl_matrix_set(X,3,1,1);

    gsl_matrix_set(X,4,0,7.3);
    gsl_matrix_set(X,4,1,2.1);

    gsl_matrix_set(X,5,0,5);
    gsl_matrix_set(X,5,1,4.1);

    int n=500; //Number of cases
    int m=2;    //Number of variables
    int k=50;

    //Number of groups
    gsl_matrix *B = gsl_matrix_alloc(n,m);
    int i,j;
    double val;
    int maxiter=100;
    double * data = (double*) malloc(n* m * sizeof(double));
    for (i=0;i<n;i++){
        for (j=0;j<m;j++){
            val= ((float)rand())/(float)RAND_MAX;
            gsl_matrix_set(B,i,j,val);
//            data[i * m + j]=val*(i+j);
        }
    }

    kmeans_data* kmd = gsl_kmeans(B,k, 100);

    // save B,centroids and print indexes
    FILE *train_data = fopen("train_data.txt","w");
    gsl_matrix_printf(B,train_data);
    fclose(train_data);

    FILE *train_center = fopen("train_center.txt","w");
    gsl_matrix_printf(kmd->centers, train_center);
    fclose(train_center);

    FILE *train_assign = fopen("train_assign.txt","w");
    for (size_t i = 0; i < kmd->index->size; ++i) {
        fprintf(train_assign, "%d", gsl_vector_int_get(kmd->index,i));
        if(i != kmd->index->size -1) {
            fprintf(train_assign," ");
        }
    }

    gsl_kmeans_free(kmd);
    return 0;
}
