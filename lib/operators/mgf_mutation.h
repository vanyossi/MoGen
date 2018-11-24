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


#ifndef MOGEN_MGF_MUTATION_H
#define MOGEN_MGF_MUTATION_H

#include "mgf_global_types.h"
#include "mogen_mop.h"

typedef enum moa_mut_types {
    MUT_PBM
} MutType;

void moa_mutation_setup(Moa *moa, MutType m_type);

void PBM_mut(Individual *indv, MutationSettings mutation, Mop_limit limits);

void PBM_real(Individual *indv, MutationSettings mutation, Mop_limit limits);

void mut_bitwise(Individual *indv, MutationSettings mutation, Mop_limit limits);


#endif //MOGEN_MGF_MUTATION_H
