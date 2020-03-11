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


#ifndef MOGEN_MOP_REPORT_H
#define MOGEN_MOP_REPORT_H

#include "mgf_global_types.h"
#include "time_tools.h"

typedef struct mop_report_stats_t {
    int gens;
    int evals;
    long t_elapsed;
} MopReportStats;

struct mgf_report_header {
    int size;
    size_t alloc_size;
    char *str;
    char *cursor;
    mbool locked;
};

typedef struct mop_report_t{
   MopReportStats current;
   MopReportStats total;
   MopReportStats record;
   unsigned int run;
   struct mgf_report_header header;
   char name[64];
   char name_var[64];
   char name_obj[64];
   void (*report)(Mop *mop);
} MopReport;


void mop_start_timer(MopReport *report);

void mop_stop_timer(MopReport *report);

void mop_restart_stats(MopReportStats *stats);

void mop_init_report(Mop *mop);

void mop_report_header(MopReport *report, const char *format, ...);

#endif //MOGEN_MOP_REPORT_H
