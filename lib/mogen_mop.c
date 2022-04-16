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

/**
 * @brief
 */

#include "mogen_mop.h"

#include <stdlib.h>
#include <stdio.h>
//#include <stdint.h>
#include <string.h>

static int mgf_count_flags(unsigned int flags, int size){
    int count = 0;
    for (int i = 0; i < size; ++i) {
        if (flags & 1) {
            count++;
        }
        flags >>= 1;
    }
    return count;
}
/**
 * @brief
 * @param mop_specs
 */
Mop *mogen_mop(char *name, MopSpecs mop_specs, size_t mem_size) {
    size_t alloc_size = sizeof(Mop);
    if (alloc_size < mem_size) {
        alloc_size = mem_size;
    }
    Mop *nmop = calloc(1, alloc_size);

    strcpy(nmop->set.name, name);

    if (CheckFlag(mop_specs, MOP_CONTIGUOUS)) { nmop->set.continuous = 1; }
    if (CheckFlag(mop_specs, MOP_RESTRICTED)) { nmop->set.restricted = 1; }

    unsigned short is_mix = 0;
    if( CheckFlag(mop_specs, MOP_REAL)) { is_mix |= MOP_REAL; }
    if( CheckFlag(mop_specs, MOP_BIN))  { is_mix |= MOP_BIN; }
    if( CheckFlag(mop_specs, MOP_INT))  { is_mix |= MOP_INT; }

    if ( mgf_count_flags(is_mix,3) > 1) {
        nmop->set.type = MOP_MIX;
    } else {
        nmop->set.type = is_mix;
    };

    if ( CheckFlag(mop_specs, MOP_DYNAMIC)) nmop->set.dynamic = 1;

    nmop->report.report = NULL;
    nmop->report.header.alloc_size = 0;
    nmop->report.header.alloc_size = 1024;
    nmop->report.header.str = calloc(nmop->report.header.alloc_size, 1); // freed on header lock
    nmop->report.header.cursor = nmop->report.header.str;

    return nmop;
}

void mop_set_params(
    Mop *mop, unsigned int real, unsigned int binary, unsigned int integer, unsigned int nobjs, unsigned int ncons) {
    mop->set.xsize = real;
    mop->set.isize = integer;
    mop->set.bsize = binary;

    mop->set.nobj = nobjs;
    mop->set.ncons = ncons;
}

void mop_set_limits_ndec(
    Mop *mop, double *xmin, double *xmax, unsigned int lim_xsize, int *imin, int *imax, unsigned int lim_isize)
{
    int j;
    if (lim_xsize) {
        mop->limits.xmin = (double*) malloc(sizeof(double) * mop->set.xsize);
        mop->limits.xmax = (double*) malloc(sizeof(double) * mop->set.xsize);

        for (size_t i = 0; i < mop->set.xsize; i++) {
            j = i % lim_xsize;
            mop->limits.xmin[i] = xmin[j];
            mop->limits.xmax[i] = xmax[j];
        }
    }
    if (lim_isize) {
        mop->limits.imin = calloc(mop->set.isize, sizeof(int));
        mop->limits.imax = calloc(mop->set.isize, sizeof(int));

        for (size_t i = 0; i < mop->set.isize; i++) {
            j = i % lim_isize;
            mop->limits.imin[i] = imin[j];
            mop->limits.imax[i] = imax[j];
        }
    }
}

Individual *mogen_mop_getindv(Mop *mop, unsigned int pos) {
    // @ TODO assert if pos is bigger than size.
    if (pos < mop->pop->size) {
        return mgf_pop_get_indv(mop->pop, pos);
    } else {
        return (void*) 0;
    }
}

mbool mop_evaluate(Mop *mop, Individual *indv) {
    mop->evaluate(mop, indv);
    mop->report.current.evals++;
    mop->report.total.evals++;
    // Check if stop condition is met
    if (moa_stop(mop->solver, MGN_STOPIF_EXEC)){
        return mfalse;
    }
    return mtrue;
}

void mop_solve(Mop *mop, int steps){
    int stop = 0;
    mop_restart_stats(&mop->report.current);

    do {
        steps--;
        mop_start_timer(&mop->report);
        // if moa returns false, stop condition is met
        if (!mgf_moa_run(mop->solver))
            stop = 1;
        mop->report.current.gens++;
        mop->report.total.gens++;

        if (moa_stop(mop->solver, MGN_STOPIF_GEN)){
            stop = 1;
        }
        mop_stop_timer(&mop->report);
        // one gen done, now report
        
        if (stop) break;
    } while (steps);

//    printf("status of run %d", stop);
//    nanosleep((const struct timespec[]){{1L, 30000000L}}, NULL);

//    mop->report.total.gens += mop->report.current.gens;
//    mop->report.total.evals += mop->report.current.evals;
}

void mop_print(Mop *mop){
    printf("Mop specs are \n"
           "Cont: %d\n"
           "Res: %d\n"
           "Type: %d\n"
           "Dyn: %d\n\n"
           , mop->set.continuous
           , mop->set.restricted
           , mop->set.type
           , mop->set.dynamic);
}
