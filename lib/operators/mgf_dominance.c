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


#include "mgf_dominance.h"

/**
 * To verify dominance between two solutions
 * Returns
 *  1: if x weak dominates B
 *  0: otherwise
 *  */
int weak_dominance(double *x, double *y, int dim)
{
    int flag = 1;
    for (int i = 0; i < dim; ++i) {
        if (x[i] > y[i]) {
            flag = 0;
        }
    }
    return flag;
}

/**
 * To verify dominance between two solutions
 * @param A:    vector 'A'
 * @param B:    vector'B'
 * @param dim:  dimension of vectors
 * @return:
 *   1: if A dominates B
 * 	-1: if B dominates A
 *   0: if they are non-dominated
 *   2:	if A = B
 */
int dominance(double *A, double *B, int dim)
{
    int j, countA = 0, countB = 0;

    for (j = 0; j < dim; j++){
        if (A[j] < B[j]) {
            countA++;

        } else if (B[j] < A[j]) {
            countB++;
        }
    }
    if (countA > 0 && countB == 0){
        return 1;
    }
    if (countB > 0 && countA == 0){
        return -1;
    }
    if (countA == 0 && countB == 0){
        return 2;
    }
    return 0;
}
