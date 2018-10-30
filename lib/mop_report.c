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
