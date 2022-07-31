//
// Created by Iván Yossi on 31/07/22.
//

#ifndef MOGEN_MGN_GNUPLOT_H
#define MOGEN_MGN_GNUPLOT_H

#include <stdbool.h>
#include "mgn_types.h"

typedef struct pmgn_plot_struct mgn_plot_data;

struct pmgn_plot_struct {
    char *filename;
    char title[32];
    char xlabel[16];
    char ylabel[16];
};

void mgn_plot_open();

void mgn_plot(mgn_pop_proto *pop, mgn_plot_data *gdata);

void mgn_plot_close();

#endif //MOGEN_MGN_GNUPLOT_H
