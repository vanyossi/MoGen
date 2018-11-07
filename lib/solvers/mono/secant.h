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


#ifndef MOGEN_SECANT_H
#define MOGEN_SECANT_H

#include "mgf_mono_solver.h"

/**
 * @brief Define secant solver method
 */
#define solve(fx, x) x[1] - ((fx(x[1]) * (x[1] - x[0])) / (fx(x[1]) - fx(x[0] + .00000001)))

Moa *moa_secant(Mop *mop, double epsilon);

void moa_secant_solver(Mop *mop, Individual *indv);


#endif //MOGEN_SECANT_H
