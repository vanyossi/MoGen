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


#ifndef MOGEN_SVM_H
#define MOGEN_SVM_H

#include "mgf_moa.h"
#include "mogen_mop.h"

typedef double (moa_svm_kernel*)(double *X, double *Theta)

typedef struct moa_svm_t {
    Moa *moa;
    moa_svm_kernel kernel;
};

moa_svm(Mop* mop, double c);

moa_svm_set_kernel(Moa *moa, moa_svm_kernel kernel);


#endif //MOGEN_SVM_H
