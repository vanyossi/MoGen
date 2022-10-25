//
// Created by Iván Yossi on 29/07/22.
//

#include "mgn_rbf.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>

#include "gsl_vector_additional.h"
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
    gsl_matrix_view train_y = gsl_matrix_submatrix(trainData,0,0,trainData->size1,3);


    // print what is inside CORRECT
//    gsl_vector *ar = gsl_vector_alloc(train_x.matrix.size2);
//    for (size_t row = 0; row < train_x.matrix.size1; ++row) {
//        gsl_matrix_get_row(ar, &train_x.matrix,row);
//        gsl_vector_fprintf(stdout, ar, "%.6f");
//    }
//    printf("===============================\n");

    // need to find weights
    // define clusters
    size_t cluster_size = 5;
    kmeans_data *km = gsl_kmeans(&train_x.matrix,cluster_size, 1000);
    cluster_data_extra *kme = gsl_kmeans_calc(km);
    cluster_data cdat = {km->centers, km->k};

//    gsl_matrix *m_variance = mgn_kmeans_cluster_var(km,kme,&train_x.matrix,true);
//    gsl_matrix_printf(m_variance,stdout);

    gsl_vector *dist_sigma = mgn_kmeans_cluster_var_dist(km,kme,&train_x.matrix,true);
    gsl_matrix *mphi = mgn_rbf_create_phi(&train_x.matrix,&cdat,dist_sigma,rbf_kernel_gauss,0);
    printf("size of sigma %zu, km %zu, phi %zu, %zu tset %zu %zu\n"
           , dist_sigma->size, km->k, train_x.matrix.size1, train_x.matrix.size1
           , mphi->size1, mphi->size2);
//    gsl_matrix_printf(mphi,stdout);

    gsl_matrix *W = mgn_rbf_new_weight(mphi,&train_y.matrix,0);

    // prefidct f from W.
    gsl_matrix *y_p = gsl_matrix_alloc(train_y.matrix.size1, train_y.matrix.size2);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,mphi,W,0,y_p);

    // using testing data
    mgn_io_datagroup test_data = {0,0};
    it_read_data("./testing_set.txt", &test_data);

    if(test_data.size == 0) {
        printf("no data in file!\n");
        return 1;
    }

    gsl_matrix *m_test_data = inData_toGSLMatrix(inGroup_getListAt(&test_data,0));
    gsl_matrix_view test_x = gsl_matrix_submatrix(m_test_data,0,0,1,m_test_data->size2);

    gsl_matrix_printf(&test_x.matrix, stdout);
    gsl_matrix *y_t = gsl_matrix_alloc(1, train_y.matrix.size2);
    gsl_matrix *mphi_t = mgn_rbf_create_phi(&test_x.matrix,&cdat,dist_sigma,rbf_kernel_gauss,0);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,mphi_t,W,0,y_t);

    gsl_matrix_printf(y_t, stdout);
//    return 0; //exit


    FILE *data = fopen("y_prediction.txt","w");
    gsl_matrix_printf(y_p,data);
    fclose(data);
//    gsl_vector_fprintf(stdout, dist_sigma,"%0.8f");

    double mse = 0;
    double mse_step = 0;
    for (size_t i = 0; i < train_y.matrix.size2; ++i) {
        gsl_vector_view v_row = gsl_matrix_column(&train_y.matrix,i);
        gsl_vector_view y_rowyp = gsl_matrix_column(y_p,i);
        mse_step = mgn_math_mse(&v_row.vector,&y_rowyp.vector);
        mse += mse_step;
    }
    mse /= (double) train_y.matrix.size2;

    printf("MSE: %.6f, %.6f\n", mse, mgn_math_mse_matrix(&train_y.matrix, y_p));

    // TODO free indata
    gsl_matrix_free(trainData);
    gsl_kmeans_free(km);
    mgn_cluster_data_extra_free(kme);
    return 0;
}
