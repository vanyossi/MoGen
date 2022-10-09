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


#ifndef MOGEN_MGF_SCALARIZATION_H
#define MOGEN_MGF_SCALARIZATION_H

//#include "mgf_global_types.h"
#include "mgf_weights.h"

/**
 * @enum scalarMethod
 * @brief Scalarization methods for MOEAD
 *
 * @var SCLM_WEI        Weighted Sum
 * @var SCLM_WEI_NORM   Normalized Weighted Sum
 * @var SCLM_WEIC       Weigthed Sum Constrained
 * @var SCLM_WEIC_NORM  Normalized Weighted Sum Contrained
 * @var SCLM_TCH        Tchebycheff
 * @var SCLM_TCH_NORM   Normalized Tchebycheff
 * @var SCLM_TCHRAY     Tchebycheff along rays
 * @var SCLM_PBI        Penalty Boundary Intersection Approach (Zhang and Li, 2007)
 * @var SCLM_PBI_NORM   Normalized Penalty Boundary Intersection Approach (Zhang and Li, 2007)
 * @var SCLM_ASF        Achievement Scalarization Function (ASF) to reach the extremes of the PF
 * @var SCLM_ASF_NORM   Normalized Achievement Scalarization Function (ASF) to reach the extremes of the PF
 * @var SCLM_PSA        Pascoletti-Serafini Approach (Pascoleti and Serafini, 1984)
 * @var SCLM_MPSA       Modified Pascoletti-Serafini Approach (Eichfelder, 2009)
 * @var SCLM_LAST       Dummy enum to find last member.
 */
enum scalarMethod {
    SCLM_WEI,               //!< Weighted Sum
    SCLM_WEI_NORM,          //!< Normalized Weighted Sum
    SCLM_WEIC,              //!< Weigthed Sum Constrained
    SCLM_WEIC_NORM,         //!< Normalized Weighted Sum Contrained
    SCLM_TCH,               //!< Tchebycheff
    SCLM_TCH_NORM,          //!< Normalized Tchebycheff
//    SCLM_TCHRAY,            //!< Tchebycheff along rays
//    SCLM_PBI,               //!< Penalty Boundary Intersection Approach (Zhang and Li, 2007)
//    SCLM_PBI_NORM,          //!< Normalized Penalty Boundary Intersection Approach (Zhang and Li, 2007)
//    SCLM_ASF,               //!< Achievement Scalarization Function (ASF) to reach the extremes of the PF
//    SCLM_ASF_NORM,          //!< Normalized Achievement Scalarization Function (ASF) to reach the extremes of the PF
//    SCLM_PSA,               //!< Pascoletti-Serafini Approach (Pascoleti and Serafini, 1984)
//    SCLM_MPSA,              //!< Modified Pascoletti-Serafini Approach (Eichfelder, 2009)
    SCLM_LAST               //!< Dummy enum to find last member.
};

typedef struct mgf_scalar_method_t {
    char scalar_name[16];
    scalarization_f func;
} ScalarMethod;

//void mgf_moa_scalarization(Individual *indv, ScalarMethod scalar_method, WeightParams wparam);

void mgf_moa_set_scalarization(struct mgf_scalar_method_t *s_m, enum scalarMethod method);

double wei(int nobjs, double *f, double *w, double *z, double *n);
double wei_norm(int nobjs, double *f, double *w, double *z, double *n);

double weicons(int nobjs, double *f, double *w, double *z, double *n);
double weicons_norm(int nobjs, double *f, double *w, double *z, double *n);

double tch(int nobjs, double *f, double *w, double *z, double *n);
double tch_norm(int nobjs, double *f, double *w, double *z, double *n);

#endif //MOGEN_MGF_SCALARIZATION_H
