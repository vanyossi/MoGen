/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_weights.h"

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>



typedef struct gsl_vector_node gsl_vector_node;
typedef struct gsl_vector_list gsl_vector_list;

struct gsl_vector_node {
    gsl_vector *v;
    gsl_vector_node *next
};

struct gsl_vector_list {
    size_t size;
    gsl_vector_node *first;
    gsl_vector_node *last;
};

gsl_vector_node* gsl_vector_node_calloc(size_t size)
{
    gsl_vector_node* out = calloc(1,sizeof(*out));
    out->v = gsl_vector_calloc(size);
    return out;
}

void gsl_vector_node_free(gsl_vector_node *node)
{
    gsl_vector_free(node->v);
    free(node);
}

typedef struct pmgn_array_num array_number;

struct pmgn_array_num {
    size_t max;
    size_t size;
    size_t *val;
};

bool array_num_sum(array_number *num)
{
    bool overflow = false;
    size_t pos = 0;
    while (pos < num->size) {
        num->val[pos]++;
        if (num->val[pos] >= num->max) {
            num->val[pos] = 0;
            pos++;
        } else{
            break;
        }
    }
    if (pos >= num->size) {
        overflow = true;
    }
    return overflow;
}

array_number* array_num_alloc(size_t size, size_t max)
{
    array_number *num = calloc(1, sizeof(*num));
    num->val = calloc(size, sizeof(num->val));
    num->max = max;
    num->size = size;

    return num;
}

void array_num_free(array_number *num)
{
    free(num->val);
    free(num);
}

void gsl_vector_simplex_construct(gsl_vector_list *vl, gsl_vector *w) {
    size_t w_val = 0;
    size_t dim = 0;

    double sum = 0;

    gsl_vector_node *node = vl->last;
    size_t node_dim = node->v->size;
    size_t w_size = w->size;
    array_number *an = array_num_alloc(node_dim, w_size);

    bool run = true;
    while (run) {
        for (size_t i = 0; i < node_dim; ++i) {
            gsl_vector_set(node->v, i, gsl_vector_get(w, an->val[i]));
        }

        sum = gsl_vector_sum(node->v);
        if (fcomp(1.0, sum, 1e-4)) {
            vl->size++;
            node->next = gsl_vector_node_calloc(node_dim);
            gsl_vector_memcpy(node->next->v, node->v);
            node = node->next;
            vl->last = node;
        }
        run = !array_num_sum(an);
    }
    // delete last
    gsl_vector_node_free(node);
    vl->size--;
    array_num_free(an);

    return;
}

gsl_matrix* mgn_weight_slattice(size_t H, size_t nf)
{
    gsl_vector_list vl = {0, 0, 0};

    gsl_vector *cval = gsl_vector_calloc(H +1);
    gsl_vector_set(cval,0,1e-6);
    for (size_t i = 1; i <= H; ++i) {
        gsl_vector_set(cval,i,i/(double)H);
    }

    vl.first = gsl_vector_node_calloc(nf);
    vl.last = vl.first;
    vl.size++;

    gsl_vector_simplex_construct(&vl,cval);

    gsl_matrix *W = gsl_matrix_calloc(vl.size, nf);
    gsl_vector_node *node = vl.first;
    gsl_vector_node *node_tmp;

    for (size_t i = 0; i < vl.size; ++i) {
        gsl_matrix_set_row(W,i,node->v);
        node_tmp = node;
        node = node->next;
        gsl_vector_node_free(node_tmp);
    }
    gsl_vector_free(cval);
    return W;
}
