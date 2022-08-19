/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_counter.h"

size_t mgn_count_add(mgn_count_ciclic cc, int value)
{
    return (cc.value + value) % cc.max;
}

void mgn_count_sum(mgn_count_ciclic *cc, int value)
{
    cc->value = (cc->value + value) % cc->max;
}
