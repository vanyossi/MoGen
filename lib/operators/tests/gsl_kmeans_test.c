/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_kmeans_new.h"
#include "mgn_random.h"

#include "gsl_vector_additional.h"

#include <string.h>


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

    int n=300; //Number of cases
    int m=2;    //Number of variables
    int k=10;

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
    kmeans_data_extra *kdat = gsl_kmeans_calc(kmd);

    for (size_t l = 0; l < kdat->size; ++l) {
        size_t size = kdat->mpos[l].size;
        printf("cluster %zu: size %zu\n",l, size);
        for (size_t i1 = 0; i1 < size; ++i1) {
            printf("%u ", kdat->mpos[l].pos[i1]);
        }
        printf("\n");
    }

    gsl_kmeans_data_extra_free(kdat);


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
