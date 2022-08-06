/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_IO_H
#define MOGEN_MGN_IO_H


#include <stdbool.h>

#include "mgn_io_datalist.h"

struct _inGroup_list;

bool it_read_data(const char *filename, struct _inGroup_list* parsedData);

#endif //MOGEN_MGN_IO_H
