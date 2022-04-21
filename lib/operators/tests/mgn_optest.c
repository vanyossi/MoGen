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
#include <math.h>

#include "mgn_vector_distance.h"
#include "mgn_types.h"
#include "mgn_random.h"
#include "mgn_moa.h"

#include "mgn_gen_operator.h"

void mgn_test(char* name, bool (*func)())
{
    printf("%s \n", name);
    if(func()){
        printf("PASS!\n");
    } else {
        printf("FAIL!\n");
    }
}

bool test_vec_ordering();
bool test_genop_sbx();
bool test_genop_pbm();

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
        int *iorder = gsl_vector_qsort(&crow.vector);
        gsl_vector_int_view irank = gsl_vector_int_view_array(iorder,vecnum);
        gsl_matrix_int_set_row(drank,i,&irank.vector);
        free(iorder);
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

    mgn_test("crossover", test_genop_sbx);

    for (int i = 0; i < 10; ++i) {
        test_genop_pbm();
    }

    gsl_matrix_free(wei);

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
//        printf("order %d ", order[i]);
        result = (orderdata[i] == vnorder.vector.data[i]);
    }
    printf("\n");
    free(order);
    return result;
}

// needs better design to deal with stochastic
bool test_genop_sbx(){
    double n = 30;
    double B;
    size_t size = 4;

    mgn_ga_sets ga_probs = {0.9, 0.01, NULL, NULL};
    ga_probs.mut_llim = calloc(size, sizeof(ga_probs.mut_llim));
    ga_probs.mut_ulim = calloc(size, sizeof(ga_probs.mut_ulim));
    for (size_t i = 0; i < size; ++i) {
        ga_probs.mut_llim[i] = 0;
        ga_probs.mut_ulim[i] = 1;
    }

    double p1[4] = {1, 1, 1, 1};
    double p2[4] = {2, 2, 2, 2};
    double c1[4] = {0, 0, 0, 0};
    double c2[4] = {0, 0, 0, 0};

    printf("%.8f %g\n", pow(-2,-0.9), mgn_pow(-2,-0.5));
    printf("%.8f %g\n", pow(2,-0.5), mgn_pow(2,-0.5));
    printf("%.16f \n", 0.5 * M_PI);
//    pow(2,-0.9) * cos(-0.9 * M_PI);
//    return true;

    for (int i = 0; i < 1; ++i) {
        mgn_genop_sbx(n,p1,p2,c1,c2,size,&ga_probs);
        B = fabs( (c2[0] - c1[0]) / (p2[0] - p1[0]) );
        for (size_t j = 0; j < size; ++j) {
            printf("%.6f %.6f %.6f %.6f\n", p1[j], p2[j], c1[j], c2[j]);
//            p1[j] = c1[j];
//            p2[j] = c2[j];
        }
        printf("\n");
    }
    rnd_getUniform();

    for (int i = 0; i < 1; ++i) {
        mgn_genop_sbx_alt(n,p1,p2,c1,c2,size,&ga_probs);
        B = fabs( (c2[0] - c1[0]) / (p2[0] - p1[0]) );
        for (size_t j = 0; j < size; ++j) {
            printf("%.6f %.6f %.6f %.6f\n", p1[j], p2[j], c1[j], c2[j]);
        }
        printf("\n");
    }

//    for (size_t i = 0; i < size; ++i) {
//        printf("c1;c2:: %.6f, %.6f\n", c1[i], c2[i]);
//    }
//    printf("Bval: %.6f\n", B);

    return (fabs(B - 1) < 0.1);
}

bool test_genop_pbm()
{
    double n = 30;
    size_t size = 1;
    double p[1] = {3};
    double lb[1] = {1};
    double ub[1] = {8};

    mgn_genop_pbm(n,0.01,p,lb,ub,size);

    printf("np: %.6f\n", p[0]);

    return true;
}


