//
// Created by Iván Yossi on 17/04/22.
//

#include "../mgn_weights.h"

#include "mgn_gnuplot.h"

void test_print_matrix(gsl_matrix *W);

void test_wei_das_dennis();

int main() {

    test_wei_das_dennis();

    return 0;
}


void test_wei_das_dennis()
{
    gsl_matrix *W = mgn_weight_slattice(60, 3);
    printf("Permutation test, size: %zu, %zu\n", W->size1, W->size2);

    mgn_plot_open();
    mgn_plot_matrix_2d(W,"wvec_comb2", "wei",0);
    mgn_plot_close();

        test_print_matrix(W);

    gsl_matrix_free(W);
}

void test_print_matrix(gsl_matrix *W)
{
    for (size_t i = 0; i < W->size1; ++i) {
        double sum = 0;
        for (size_t j = 0; j < W->size2; ++j) {
            sum += gsl_matrix_get(W,i,j);
            printf("%.6f ", gsl_matrix_get(W,i,j));
        }
        printf("%.6f\n", sum);
    }
}
