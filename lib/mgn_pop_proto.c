/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_pop_proto.h"


void mgn_pop_print(mgn_pop_proto *popp, FILE *stream)
{
    gsl_vector *x;
    gsl_vector *f;
    gsl_vector *g;

    char *format = calloc(64, sizeof(*format));
    fprintf(stream,"# ");

    for (size_t fi = 0; fi < popp->iparams.f_size; ++fi) {
        asprintf(&format, "<obj%zu> ", fi);
        fprintf(stream,"%s",format);
    }

    for (size_t gi = 0; gi < popp->iparams.g_size; ++gi) {
        asprintf(&format, "<con#%zu> ", gi);
        fprintf(stream,"%s",format);
    }

    for (size_t xi = 0; xi < popp->iparams.x_size; ++xi) {
        asprintf(&format, "<var#%zu>", xi);
        fprintf(stream,"%s",format);
        if (xi != popp->iparams.x_size -1) {
            fprintf(stream," ");
        }
    }
    fprintf(stream,"\n");

    for (size_t i = 0; i < popp->size; ++i) {
        x = popp->ops->get_iparams(popp->get(popp,i)).x;
        f = popp->ops->get_iparams(popp->get(popp,i)).f;
        g = popp->ops->get_iparams(popp->get(popp,i)).g;

        for (size_t fi = 0; fi < f->size; ++fi) {
            fprintf(stream,"%e ",f->data[fi]);
        }

        for (size_t gi = 0; gi < g->size; ++gi) {
            fprintf(stream,"%e ",g->data[gi]);
        }

        for (size_t xi = 0; xi < x->size; ++xi) {
            fprintf(stream,"%e",x->data[xi]);
            if (xi != x->size -1) {
                fprintf(stream," ");
            }
        }
        fprintf(stream,"\n");
    }
}
