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


#ifndef MOGEN_TIME_TOOLS_H
#define MOGEN_TIME_TOOLS_H

#include <time.h>

#if defined(__MACH__)
#include <sys/time.h>
#endif

#if !defined(CLOCK_REALTIME)
#define CLOCK_REALTIME 0
#endif

void mogen_time_to(struct timespec *ts);

long mogen_cac_elapsed_time_ts(struct timespec *start, struct timespec *end);

long mogen_cac_elapsed_time_l(long start, long end);

void mogen_add_timeto(struct timespec *t1, struct timespec *t2);

long mogen_timespec2long(struct timespec *t);

double mogen_time_ms2sec(long time);

#endif //MOGEN_TIME_TOOLS_H
