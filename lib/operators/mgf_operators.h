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


#ifndef MOGEN_MGF_OPERATORS_H
#define MOGEN_MGF_OPERATORS_H

#include "mogen_mop.h"

typedef void (*mgf_op_crossover)(Mop*, Individual*, Individual*, Individual*, Individual*);

typedef struct mgf_operators {
    mgf_op_crossover crossover;
} Operators;

void mgf_opset_crossover(mgf_op_crossover crossover);

#endif //MOGEN_MGF_OPERATORS_H
