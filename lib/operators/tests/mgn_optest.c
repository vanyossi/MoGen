//
// Created by Iván Yossi on 15/04/22.
//

//#ifndef MOGEN_MGN_OPTEST_H
//#define MOGEN_MGN_OPTEST_H
//
//#endif //MOGEN_MGN_OPTEST_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "mgn_vector_distance.h"
#include "mgn_types.h"
#include "mgn_random.h"

void mgn_test(char* name, bool (*func)())
{
    printf("%s ", name);
    if(func){
        printf("PASS!\n");
    } else {
        printf("FAIL!\n");
    }
}

bool test_vec_ordering();

int main(int argc, char const *argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    printf("test ran!\n");
    int vecnum = 4;
    int vecdim = 2;

    gsl_matrix *wei = gsl_matrix_alloc(vecnum,vecdim);
    for (int i = 0; i < vecnum; ++i) {
        gsl_vector *temp_vec = gsl_vector_alloc(vecdim);
        gsl_matrix_get_row(temp_vec,wei, i);
//        gsl_vector_set_all(temp_vec,i);
        gsl_vector_set_seq(temp_vec);
        gsl_vector_add_constant(temp_vec, i);
        gsl_vector_scale(temp_vec, i);
        gsl_vector_fprintf(stdout,temp_vec,"%0.8f");
        printf("=====\n");
        gsl_matrix_set_row(wei,i,temp_vec);
        gsl_vector_free(temp_vec);
    }
//    gsl_matrix_fprintf(stdout, wei, "%g");
    printf("sizes %zu, %zu\n", wei->size1, wei->size2);
    gsl_matrix *dist = gsl_vector_distance_matrix(wei,2);
    gsl_matrix_int *drank = gsl_matrix_int_alloc(dist->size1, dist->size2);
    for (int i = 0; i < vecnum; ++i) {
        gsl_vector_view crow = gsl_matrix_column(dist,i);
        gsl_vector_int_view irank = gsl_vector_int_view_array(gsl_vector_qsort(&crow.vector),vecnum);
        gsl_matrix_int_set_row(drank,i,&irank.vector);
    }
    gsl_matrix_int_fprintf(stdout, drank, "%d");
    gsl_matrix_int_free(drank);
    gsl_matrix_free(dist);

    // ordering test
    mgn_test("ordering", test_vec_ordering);

    double data[6] = {5, 4, 3, 2, 1, 100};
    gsl_vector_view vdat = gsl_vector_view_array(data, 6);
    int *order = gsl_vector_qsort(&vdat.vector);
    for (size_t i = 0; i < 6; ++i) {
        printf("data; %d %.5f\n", order[i], data[i]);
    }
    free(order);



//    gsl_vector *Vec = calloc(vecnum, sizeof(gsl_vector));
//    for (int i = 0; i < vecnum; ++i) {
//        Vec[i] = *gsl_vector_alloc(vecdim);
//        gsl_vector_set_all(&Vec[i],i);
//    }
////    for (int i = 0; i < vecnum; ++i) {
////        gsl_vector_fprintf(stdout,&Vec[i],"%.5f");
////    }
//    printf("====\n");
//    gsl_matrix *dist2 = gsl_vector_distance_matrix1(Vec,vecnum,3);
////    gsl_matrix_fprintf(stdout, dist2, "%.6f");
//    gsl_matrix_free(dist2);


    return 0;
}

bool test_vec_ordering() {
    bool result = true;
    size_t vecdim = 6;
    double testdata[6] = {4,3, 2,5, 1, 6};
    double orderdata[6] = {1,2,3,4,5,6};

    gsl_vector_view vnorder = gsl_vector_view_array(testdata,vecdim);
    int *order = gsl_vector_qsort(&vnorder.vector);

    for (size_t i = 0; i < vecdim; ++i) {
        printf("order %d ", order[i]);
        result = (orderdata[i] == vnorder.vector.data[i]);
    }
    printf("\n");
    free(order);
    return result;
}
