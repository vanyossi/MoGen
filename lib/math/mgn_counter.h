/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_COUNTER_H
#define MOGEN_MGN_COUNTER_H

#include <stdlib.h>

typedef struct mgnp_ccounter mgn_count_ciclic;

struct mgnp_ccounter {
    size_t max;
    size_t value;
};

size_t mgn_count_add(mgn_count_ciclic cc, int value);

void mgn_count_sum(mgn_count_ciclic *cc, int value);

#endif //MOGEN_MGN_COUNTER_H
