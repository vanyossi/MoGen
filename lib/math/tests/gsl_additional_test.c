/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "gsl_vector_additional.h"

int main(int argc, char const *argv[]) {
    gsl_vector *iter = gsl_vector_alloc(6);
    for (size_t i = 0; i < iter->size; ++i) {
        gsl_vector_set(iter,i,i);
    }

    gsl_matrix *rep = gsl_matrix_alloc(3,iter->size);
    gsl_vector_repeat(iter, 3, rep);
    gsl_vector_view vec = gsl_matrix_row(rep,1);
    gsl_vector_fprintf(stdout, &vec.vector, "%.6f");

    gsl_matrix_free(rep);
    gsl_vector_free(iter);
    return 0;
}
