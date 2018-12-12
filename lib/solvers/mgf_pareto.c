/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "mgf_pareto.h"

#include <float.h>
#include <math.h>

#include "mgf_population.h"

void mgf_pareto_update_z(Individual *indv, double *z, int dim)
{
    for (int i = 0; i < dim; ++i) {
        z[i] = fmin(indv->f[i], z[i]);
    }
}

/**
 * Calculate the Nadir value of population
 * @param pop MoeazPop population
 * @param nadir output array of nadir values
 * @param nobj number of objectives
 */
void mgf_pop_get_nadir(MoeazPop *pop, double *nadir, int nobj)
{
    unsigned int i, j;

    for (i = 0; i < nobj; ++i) {
        nadir[i] = -DBL_MAX;
    }

    for (i = 0; i < pop->size; ++i) {
        for (j = 0; j < nobj; ++j) {
            nadir[j] = fmax(mgf_pop_get_indv(pop,i)->f[j], nadir[j]);
        }
    }
}
