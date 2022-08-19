#ifndef _MGN_MOEADRBF_SOLVER_
#define _MGN_MOEADRBF_SOLVER_


#include "mgn_types.h"
#include "mgn_moa.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


typedef struct mgn_features_moeadrbf mgn_moeadrbff;


#define pmgn_moeadrbf_templatep() \
    size_t max_eval;              \
    size_t nt; /* size of tset population */ \
    gsl_matrix *m_w; /* matrix of weight vectors */ \
    mgn_popl *solution;


struct mgn_features_moeadrbf {
    pmgn_moeadrbf_templatep()
};


// public functions and structures
// moeadrbf_data is for public use
// used for initialization data retrieval and else
mgnMoa* mgn_moa_moeadrbf_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    );

void mgn_moa_moeadrbf_free(mgnMoa* moeadrbf);

//void mgn_de_setmop(
//    mgnMoa *mrbf
//    , mgnMop *mop
//    , void* todo);


void mgn_moa_moeadrbf_free(mgnMoa *mrbf);

void mgn_moa_moeadrbf_init(mgnMoa* moeadrbf);

#endif // _MGN_MOEADRBF_SOLVER_
