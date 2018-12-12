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


#ifndef MOGEN_ZDT_H
#define MOGEN_ZDT_H

#include "mgf_global_types.h"
#include "mogen_mop.h"

typedef enum mop_zdt_e {
    ZDT1,
    ZDT2,
    ZDT3,
    ZDT4,
    ZDT6,
    ZDTM1
} ZDTVariant;

Mop *mop_zdt(ZDTVariant zdt);

#endif //MOGEN_ZDT_H
