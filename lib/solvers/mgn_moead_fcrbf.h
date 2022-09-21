#ifndef _MGN_MOEAD_FCRBF_SOLVER_
#define _MGN_MOEAD_FCRBF_SOLVER_


#include "mgn_types.h"
#include "mgn_moa.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


typedef struct mgn_features_moead_fcrbf mgn_moead_fcrbff;


#define pmgn_moead_fcrbf_templatep() \
    size_t max_eval;              \
    size_t nt; /* size of tset population */ \
    gsl_matrix *m_w; /* matrix of weight vectors */ \
    mgn_popl *solution;


struct mgn_features_moead_fcrbf {
    pmgn_moead_fcrbf_templatep()
};


// public functions and structures
// moead_fcrbf_data is for public use
// used for initialization data retrieval and else
mgnMoa* mgn_moa_moead_fcrbf_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    , size_t Nw /* size of internal weight vector */
);

void mgn_moa_moead_fcrbf_init(mgnMoa* moead_fcrbf);

void mgn_moa_moead_fcrbf_free(mgnMoa* moead_fcrbf);


//void mgn_moead_fcrbf_mdl_set_psize(mgnMoa* moead_fcrbf);



#endif // _MGN_MOEAD_FCRBF_SOLVER_
