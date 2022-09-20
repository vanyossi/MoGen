//
// Created by Iván Yossi on 31/07/22.
//

#ifndef MOGEN_MGN_GNUPLOT_H
#define MOGEN_MGN_GNUPLOT_H

#include <stdbool.h>
#include <gsl/gsl_matrix.h>

#include "mgn_types.h"

typedef struct pmgn_plot_struct mgn_plot_data;

struct pmgn_plot_struct {
    char *filename;
    char *title;
    char xlabel[16];
    char ylabel[16];
    float lxrange;
    float uxrange;
    float lyrange;
    float uyrange;
};

void mgn_plot_open();

void mgn_plot(mgn_pop_proto *pop, mgn_plot_data *gdata);

void mgn_plot_fast(mgn_pop_proto *i_pop, char* filename, char* title);

void mgn_plot_matrix_2d(gsl_matrix *m_mat
                        , const char* filename
                        , const char* title
                        , size_t *idx);

void mgn_plot_close();

#endif //MOGEN_MGN_GNUPLOT_H
