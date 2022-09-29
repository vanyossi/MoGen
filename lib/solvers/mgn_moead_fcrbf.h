#ifndef MGN_MOEAD_FCRBF_SOLVER
#define MGN_MOEAD_FCRBF_SOLVER

#include "mgn_moead-rbf_common.h"


// public functions and structures
// moeadrbf_data is for public use
// used for initialization data retrieval and else
mgnMoa* mgn_moa_moeadrbf_fc_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    , size_t Nw /* size of internal weight vector */
);

void mgn_moa_moeadrbf_fc_init(mgnMoa* moeadrbf);

void mgn_moa_moeadrbf_fc_free(mgnMoa* moeadrbf);


#endif // MGN_MOEAD_FCRBF_SOLVER
