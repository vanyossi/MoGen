//
// Created by Iván Yossi on 17/04/22.
//

#include "../mgn_weights.h"

void test_print_matrix(gsl_matrix *W);

void test_wei_combination();
void test_wei_permutation();

int main() {
//    gsl_matrix *W = mgn_weight_slattice_comb(6, 2);
//    gsl_matrix_fprintf(stdout,W,"%.6f");
//    gsl_matrix_free(W);

    test_wei_combination();
    test_wei_permutation();

//    for (size_t i = 0; i < 30; ++i) {
//        printf("%zu -> %zu\n", i, (i+1) % 5);
//    }
    return 0;
}

void test_wei_combination()
{
    gsl_matrix *W = mgn_weight_slattice_comb(3, 2);
    puts("Combination test");
    test_print_matrix(W);

    gsl_matrix_free(W);
}

void test_wei_permutation()
{
    gsl_matrix *W = mgn_weight_slattice_perm(3, 2);
    puts("Permutation test");
    test_print_matrix(W);
//    for (size_t i = 0; i < W->size1; ++i) {
//        for (size_t j = 0; j < W->size2; ++j) {
//            printf("%.6f ", gsl_matrix_get(W,i,j));
//        }
//        printf("\n");
//    }
    gsl_matrix_free(W);
}

void test_print_matrix(gsl_matrix *W)
{
    for (size_t i = 0; i < W->size1; ++i) {
        for (size_t j = 0; j < W->size2; ++j) {
            printf("%.6f ", gsl_matrix_get(W,i,j));
        }
        printf("\n");
    }
}
