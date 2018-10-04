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
#include <string.h>

#define CheckSpec(flag) (mop_specs & flag) == flag

/**
 * @brief
 * @param mop_specs
 */
Mop *mogen_mop(char *name, MopSpecs mop_specs) {

    Mop *nmop = calloc(1, sizeof(Mop));
    strcpy(nmop->set.name, name);

    if (CheckSpec(MOP_CONTIGUOUS)) {
        nmop->set.continuous = 1;
    }
    if (CheckSpec(MOP_RESTRICTED)) {
        nmop->set.restricted = 1;
    }

    unsigned short is_mix = 0;
    if( CheckSpec(MOP_REAL)) {
        is_mix |= MOP_REAL;
    }
    if( CheckSpec(MOP_BIN)) {
        is_mix |= MOP_BIN;
    }
    nmop->set.type = is_mix;

    if ( CheckSpec(MOP_DYNAMIC)) {
        nmop->set.dynamic = 1;
    }

    return nmop;
}

MoeazIndv *mogen_mop_getIndv(Mop *mop, unsigned int pos) {
    return &mop->pop->indv[pos];
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
