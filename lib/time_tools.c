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

#include "time_tools.h"

void mogen_time_to(struct timespec *ts){
#if CLOCK_REALTIME
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_MONOTONIC, ts);
#endif
}

long mogen_cac_elapsed_time(struct timespec *start, struct timespec *end){
    long ms_diff;
    ms_diff = (end->tv_nsec - start->tv_nsec)/1000; //microsegundos
    ms_diff += (end->tv_sec - start->tv_sec) * 1000000;
    return ms_diff;
//    double deltat_s  = end->tv_sec - start->tv_sec;
//    double deltat_ns = end->tv_nsec - start->tv_nsec;
//    return deltat_s + deltat_ns*1e-9;
}

long mogen_cac_elapsed_time_l(long start, long end){
    return end - start;
}

void mogen_add_timeto(struct timespec *t1, struct timespec *t2){
    t1->tv_sec += t2->tv_sec;
    t1->tv_nsec += t2->tv_nsec;
}

long mogen_timespec2long(struct timespec *t) {
    long ms_time;
    ms_time = (long)(t->tv_sec * 1e6);
    ms_time += (long)(t->tv_nsec * 1e-3);
    return ms_time;
}

double mogen_time_ms2sec(long time){
    return time * 1e-6;
}
