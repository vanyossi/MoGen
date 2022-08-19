//
// Created by Iván Yossi on 31/07/22.
//

#include "mgn_gnuplot.h"
#include "gnuplotbin.h"

#include <stdio.h>

#include "population.h"

static FILE *gnuplot;

void mgn_plot_open()
{
    gnuplot = popen(GNUPLOT, "w");
}

void mgnp_plot_w_headers(mgn_plot_data* gdata)
{
    fprintf(gnuplot, "\nreset\n");
    fprintf(gnuplot, "set terminal png\n");
    fprintf(gnuplot, "set output \"%s.png\"\n", gdata->filename);
    fprintf(gnuplot, "plot '-' title \"%s\"\n", gdata->title);
}

void mgnp_plot_w_closure()
{
    fprintf(gnuplot, "e\n");
//    fprintf(gnuplot, "set output\n");
    fflush(gnuplot);
}

// TODO moeaz uses a report_file
// fixed to 2D
void mgn_plot(mgn_pop_proto *popp, mgn_plot_data *gdata)
{
    mgn_pop* pop = (mgn_pop*)popp;
    fprintf(gnuplot, "\nreset\n");
    fprintf(gnuplot, "set terminal png\n");
    fprintf(gnuplot, "set output \"%s.png\"\n", gdata->filename);
//    fprintf(gnuplot, "set xrange[%f:%f]\n", gdata->lxrange, gdata->uxrange);
//    fprintf(gnuplot, "set yrange[%f:%f]\n", gdata->lyrange, gdata->uyrange);
    fprintf(gnuplot, "plot '-' title \"non dominated\"\n");

    for (size_t i = 0; i < pop->size; i++) {
        double *f = pop->ops->get_iparams(pop->get(pop,i)).f->data;
        fprintf(gnuplot, "%g %g\n", f[0], f[1]);
    }
    fprintf(gnuplot, "e\n");
//    fprintf(gnuplot, "set output\n");
    fflush(gnuplot);
}

void mgn_plot_fast(mgn_pop_proto *i_pop, char* filename, char* title)
{
    mgn_pop* pop = (mgn_pop*)i_pop;
    mgn_plot_data gdata = {filename, title
                           , "", "", 0,0,0,0};
    mgnp_plot_w_headers(&gdata);

    for (size_t i = 0; i < pop->size; i++) {
    double *f = pop->ops->get_iparams(pop->get(pop,i)).f->data;
        fprintf(gnuplot, "%g %g\n", f[0], f[1]);
    }

    mgnp_plot_w_closure();
}

void mgn_plot_matrix_2d(gsl_matrix *m_mat, char* filename, char* title, size_t *idx)
{
    if (idx == 0) {
        size_t idx_l[2] = {0, 1};
        idx = idx_l;
    }
    mgn_plot_data gdata = {filename, title
        , "", "", 0,0,0,0};
    mgnp_plot_w_headers(&gdata);

    for (size_t i = 0; i < m_mat->size1; ++i) {
        gsl_vector_view vv = gsl_matrix_row(m_mat,i);
        fprintf(gnuplot, "%g %g\n", vv.vector.data[idx[0]], vv.vector.data[idx[1]]);
    }

    mgnp_plot_w_closure();
}

void mgn_plot_close()
{
    fclose(gnuplot);
    gnuplot = NULL;
}
