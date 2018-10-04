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
 * @file multiarray.h
 * @author Iván Santa María
 * @date 11 Sept 2018
 * @brief Multiarray, double, int
 *
 * Long description
 *
 */

#ifndef MULTIARRAY_H
#define MULTIARRAY_H

#include <stdlib.h>


typedef struct {
    int         size;       //!< Multiarray members size
    int         muatype;    //!< Defines type of multiarray (real, bin, mixed)
    int*        type_idx;   //!< Type index identify which position are binary
    double*     real;       //!< Double array for real values
    int*        bin;        //!< Int array, each bit represents one position
} Multiarray;


/**
 * @brief Create new Multiarray structure
 * @param ma
 * @param type
 */
void mua_multiarray(Multiarray* ma, int type);

/**
 * @brief Assign double value to position
 * @param ma
 * @param pos
 * @param value
 */
void mua_set_double(Multiarray* ma, unsigned int pos, double value);


/**
 * @brief Assing binary value to position
 * @param ma
 * @param pos
 * @param value
 */
void mua_set_int(Multiarray* ma, unsigned int pos, int value);


/**
 * @brief Return value from positon pos
 * @param ma
 * @param pos
 * @return
 */
void* mua_value_at(Multiarray* ma, unsigned int pos);


/**
 * @brief Return true if value at position is binary
 * @param ma
 * @param pos
 * @return
 */
int mua_value_isbin(Multiarray* ma, unsigned int pos);

/**
 * @brief Free memory allocated by MultiArray struct
 * @param ma
 */
void mua_multiarray_del(Multiarray *ma);

#endif /* MULTIARRAY_H */
