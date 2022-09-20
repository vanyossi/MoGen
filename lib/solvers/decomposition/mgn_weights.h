//
// Created by Iván Yossi on 17/04/22.
//

#ifndef MOGEN_MNG_WEIGHTS_H
#define MOGEN_MNG_WEIGHTS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#define fcomp(a,b,ep) fabs(a - b) < ep


gsl_matrix* mgn_weight_slattice(size_t H, size_t nf);


#endif //MOGEN_MNG_WEIGHTS_H
