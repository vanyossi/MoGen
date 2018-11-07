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


#ifndef MOGEN_GLOBAL_TYPES_H
#define MOGEN_GLOBAL_TYPES_H

typedef enum mop_specs_e MopSpecs;

typedef struct mop_t Mop;
typedef struct mop_base_t Mop_base;
typedef struct mop_extra_t Mop_extra;
typedef struct mop_limit_t Mop_limit;


typedef struct mogen_moa_t Moa;
typedef enum moa_stop_criterion_t MoaStopCriterion;
typedef enum mogen_moa_types_e MoaTypes;

typedef struct indv_type_t IndvidualType;

typedef struct indv_t Individual;

typedef struct mgf_pop_t MoeazPop;


// Useful defines
#define UNUSED(expr) do { (void)(expr); } while (0)

#endif //MOGEN_GLOBAL_TYPES_H
