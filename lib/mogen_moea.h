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


#ifndef MOGEN_MOGEN_MOEA_H
#define MOGEN_MOGEN_MOEA_H


union mogen_moa_type_t;
struct mop_t;

typedef struct moea_t {
    double cross_index; //!< Crossover index
    double cross_prob;  //!< Crossover probability of happening
    double mut_index;   //!< Mutation index
    double mut_prob;    //!< Mutation probability of happening
} Moea;


struct mogen_moa_t *moa_moea(struct mop_t *mop, union mogen_moa_type_t moea, char *name);

#endif //MOGEN_MOGEN_MOEA_H


