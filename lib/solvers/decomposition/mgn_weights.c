/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_weights.h"

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>

#include "mgn_random.h"

// campelo2018moeadr
gsl_matrix* mgn_weight_slattice_comb(size_t H, size_t nf)
{
    gsl_combination *comb = gsl_combination_calloc(H,nf);

    size_t rows = (int)round(gsl_sf_choose(H,nf));

    gsl_vector *cval = gsl_vector_alloc(H);
    for (size_t i = 1; i <= H; ++i) {
        gsl_vector_set(cval,i-1,i/(double)H);
    }

//    printf("comb size: %zu  %zu", comb->k, comb->n);

    gsl_matrix *W = gsl_matrix_alloc(rows,nf);

//    gsl_vector_ulong *vcc = gsl_vector_ulong_alloc(nf);
    for (size_t i = 0; i < rows; ++i) {
//        gsl_vector_ulong_view vc = gsl_vector_ulong_view_array(comb->data, nf);
//        gsl_vector_ulong_memcpy(vcc, &vc.vector);
//
//        gsl_ran_shuffle(rnd_get_generator(),vcc->data, nf, sizeof(size_t));

        for (size_t j = 0; j < nf; ++j) {
            gsl_matrix_set(W,i,j,gsl_vector_get(
                cval, gsl_combination_get(comb,j)));
        }
        gsl_combination_next(comb);
    }
//    gsl_vector_ulong_free(vcc);
    gsl_combination_free(comb);
    gsl_vector_free(cval);
    return W;
}

gsl_matrix* mgn_weight_slattice_perm(size_t H, size_t nf)
{
    size_t h_with_zero = H + 1;
    gsl_combination *comb = gsl_combination_calloc(h_with_zero,nf);
    gsl_permutation *perm = gsl_permutation_calloc(nf);
    size_t rows = (int)round(gsl_sf_choose(h_with_zero,nf)) * perm->size;

    gsl_vector *cval = gsl_vector_alloc(h_with_zero);
    gsl_vector_set(cval,0,1e-6);
    for (size_t i = 1; i <= H; ++i) {
        gsl_vector_set(cval,i,i/(double)H);
    }

//    printf("comb size: %zu  %zu", comb->k, comb->n);
    int feasable_c = 0;
    int *feasable = calloc(rows,sizeof(*feasable));
    gsl_matrix *W = gsl_matrix_alloc(rows,nf);

    gsl_vector_ulong *pcomb = gsl_vector_ulong_alloc(nf);
    for (size_t i = 0; i < rows; ++i) {
        gsl_vector_view crow = gsl_matrix_row(W,i);
        gsl_vector_ulong_view combv = gsl_vector_ulong_view_array(comb->data, nf);
        for (size_t j = 0; j < perm->size; ++j) {
            gsl_vector_ulong_set(pcomb,j, combv.vector.data[gsl_permutation_get(perm,j)]);
        }
        for (size_t j = 0; j < nf; ++j) {
            gsl_vector_set(&crow.vector,j,gsl_vector_get(
                cval, gsl_vector_ulong_get(pcomb,j)));
//            gsl_matrix_set(W,i,j,gsl_vector_get(
//                cval, gsl_vector_ulong_get(pcomb,j)));
        }
        if ( fcomp(1.0,gsl_vector_sum(&crow.vector),1e-4) ) {
            feasable[i] = (int)i;
            feasable_c++;
        } else {
            feasable[i] = -1;
        }

        gsl_permutation_next(perm);

        if ((i+1) % perm->size == 0) {
            gsl_permutation_init(perm);
            gsl_combination_next(comb);
        }
    }

//    for (size_t i = 0; i < W->size1; ++i) {
//        gsl_vector_view row = gsl_matrix_row(W,i);
//        printf("s %.6f\n", gsl_vector_sum(&row.vector));
//    }

    // select only w vectors that sum 1
    gsl_matrix *Wtrim = gsl_matrix_alloc(feasable_c, nf);
//    printf("feas %d\n\n", feasable_c);
    size_t j = 0;
    for (size_t i = 0; i < rows; ++i) {
        if(feasable[i] > -1) {
            gsl_vector_view cr = gsl_matrix_row(W,feasable[i]);
            gsl_matrix_set_row(Wtrim,j, &cr.vector);

            j++;
        }
    }

    free(feasable);
    gsl_vector_ulong_free(pcomb);
    gsl_permutation_free(perm);
    gsl_combination_free(comb);
    gsl_vector_free(cval);
    gsl_matrix_free(W);
    return Wtrim;
}
