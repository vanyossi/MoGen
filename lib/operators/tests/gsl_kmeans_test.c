/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_kmeans_new.h"
#include "mgn_random.h"

#include "gsl_vector_additional.h"
#include "mgn_io.h"

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

    /**
     * @result
        66.0     55.0
        48.0     52.0
        27.75    55.0
        56.75    35.75
        43.2     16.7
        30.8333  74.6667
     */
    char* fname = "kmeans_test.txt";
    struct _inGroup_list io_data_fp = {0, 0};
    it_read_data(fname,&io_data_fp); // alters size of pl_a
    gsl_matrix *B = inData_toGSLMatrix(inGroup_getListAt(&io_data_fp,0));

    k = 3;

    kmeans_data* kmd = gsl_kmeans(B,k, 100);
    cluster_data_extra *kdat = gsl_kmeans_calc(kmd);
    printf("max iter: %lu\n", kmd->iter);

    for (size_t l = 0; l < kdat->size; ++l) {
        size_t size = kdat->mpos[l].size;
        printf("cluster %zu: size %zu\n",l, size);
        for (size_t i1 = 0; i1 < size; ++i1) {
            printf("%u ", kdat->mpos[l].pos[i1]);
        }
        printf("\n");
    }


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
    fclose(train_assign);

    mgn_cluster_data_extra_free(kdat);
    gsl_kmeans_free(kmd);
    gsl_matrix_free(B);
    gsl_matrix_free(X);
    rnd_gen_free();
    return 0;
}
