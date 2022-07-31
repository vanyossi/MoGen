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

// TODO moeaz uses a report_file
// fixed to 2D
void mgn_plot(mgn_pop_proto *popp, mgn_plot_data *gdata)
{
    mgn_pop* pop = (mgn_pop*)popp;
    fprintf(gnuplot, "\nreset\n");
    fprintf(gnuplot, "set terminal png\n");
    fprintf(gnuplot, "set output \"%s.png\"\n", gdata->filename);
    fprintf(gnuplot, "plot '-' title \"non dominated\"\n");

    for (size_t i = 0; i < pop->size; i++) {
        double *f = pop->ops->get_iparams(pop->get(pop,i)).f->data;
        fprintf(gnuplot, "%g %g\n", f[0], f[1]);
    }
    fprintf(gnuplot, "e\n");
//    fprintf(gnuplot, "set output\n");
    fflush(gnuplot);
}

void mgn_plot_close()
{
    fclose(gnuplot);
}
