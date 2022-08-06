//
// Created by Iván Yossi on 29/07/22.
//

#include "mgn_rbf.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "mgn_io.h"

int main(int argc, char const *argv[]) {

    mgn_io_datagroup inputData = {0,0};
    it_read_data("./training_set.txt", &inputData);

    if(inputData.size == 0) {
        printf("no data in file!\n");
        return 1;
    }

    // generate Train data
    gsl_matrix *trainData = inData_toGSLMatrix(inGroup_getListAt(&inputData,0));
    gsl_matrix_view train_x = gsl_matrix_submatrix(trainData,0,3,trainData->size1,trainData->size2-3);

    // print what is inside CORRECT
//    gsl_vector *ar = gsl_vector_alloc(train_x.matrix.size2);
//    for (size_t row = 0; row < train_x.matrix.size1; ++row) {
//        gsl_matrix_get_row(ar, &train_x.matrix,row);
//        gsl_vector_fprintf(stdout, ar, "%.6f");
//    }
//    printf("===============================\n");

    // define cluster sizes
    size_t cluster_size = 10;
    kmeans_data* km = gsl_kmeans(&train_x.matrix,cluster_size, 1000);

    // TODO free indata
    gsl_matrix_free(trainData);

    gsl_kmeans_free(km);
    return 0;
}
