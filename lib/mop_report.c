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


#include "mop_report.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "mogen_mop.h"

void mop_start_timer(MopReport *report){
    struct timespec now;
    mogen_time_to(&now);
    report->current.t_elapsed = mogen_timespec2long(&now);
}

void mop_stop_timer(MopReport *report){
    struct timespec stop;
    mogen_time_to(&stop);

    report->current.t_elapsed = mogen_cac_elapsed_time_l(report->current.t_elapsed, mogen_timespec2long(&stop));
    report->total.t_elapsed += report->current.t_elapsed;
}

void mop_restart_stats(MopReportStats *stats){
    stats->gens = 0;
    stats->evals = 0;
}

void mop_report_start(Mop *mop){
    if (mop == NULL || mop->solver == NULL) {
        printf("Mop or Moa not defined, cannot initialize report");
        return;
    }

    sprintf(mop->report.name_obj, "%s_%s_%u.obj", mop->set.name, mop->solver->name, mop->report.run);
    sprintf(mop->report.name_var, "%s_%s_%u.var", mop->set.name, mop->solver->name, mop->report.run);

    FILE *obj_f = fopen(mop->report.name_obj,"w");
    fwrite(mop->report.header.str, sizeof(char), (size_t)mop->report.header.size, obj_f);
    fclose(obj_f);
//    fwrite("MOGEN1", 6, 1, rep_file);

    FILE *var_f = fopen(mop->report.name_obj,"w");
    fwrite(mop->report.header.str, sizeof(char), (size_t)mop->report.header.size, var_f);
    fclose(var_f);

    mop->report.header.locked = mtrue;
    free(mop->report.header.str);
}

void mgf_update_cursor(int *size, char **str)
{
    int counter = 0;
    while ((*str)[counter] != 0) {
        counter++;
    }
    *str = (*str) + counter;
    *size += counter;
}

int mgf_realloc_buffer(size_t size, size_t *alloc_size, void **buffer)
{
    int bool_realloc = 0;
    void *new_buffer;
    size_t new_alloc_size = *alloc_size * 2;
    if (size > *alloc_size){
        new_alloc_size *= 2;
        new_buffer = realloc(*buffer, new_alloc_size);
        memcpy(new_buffer, *buffer, *alloc_size);
        free(*buffer); // delete old buffer

        *buffer = new_buffer;
        bool_realloc = 1;
    }
    *alloc_size = new_alloc_size;
    return bool_realloc;
}

void mop_report_header(MopReport *report, const char *format, ...)
{
    if (report->header.locked) {
        printf("WARN: Header closed! No new data added to file headers.");
        return;
    }
    va_list arg_list;
    va_start(arg_list, format);
    vsprintf(report->header.cursor, format, arg_list);
    mgf_update_cursor(&report->header.size, &report->header.cursor);
    va_end(arg_list);
}
