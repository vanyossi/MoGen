/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_MOEAD_RBF_COMMON_H
#define MOGEN_MGN_MOEAD_RBF_COMMON_H

#include "mgn_types.h"
#include "mgn_moead_types.h"
#include "mgn_moa.h"

#include "mgn_counter.h"
#include "mgn_initializer.h"
#include "mgn_pop_matrix.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


typedef struct mgnp_model_rbfdata mgn_model_rbfdata;


#define pmgn_moeadrbf_templatep() \
    size_t max_eval;              \
    size_t nt; /* size of tset population */ \
    gsl_matrix *m_w; /* matrix of weight vectors */ \
    mgn_popl *solution;


struct mgn_features_moeadrbf {
    pmgn_moeadrbf_templatep()
};

/*
 * Functions to overload
 */
void mgn_moa_moeadrbf_common_init(mgnMoa* moeadrbf);

void mgn_moa_moeadrbf_common_free(mgnMoa* moeadrbf);

// public functions and structures
// moeadrbf_data is for public use
// used for initialization data retrieval and else
mgnMoa* mgn_moa_moeadrbf_common_alloc(
    size_t max_eval
    , size_t nt
    , size_t k
    , gsl_matrix *W
    , mgn_popl *A
    , mgnLimit *limits
    , size_t Nw /* size of internal weight vector */
);

/*
 * Private Structures
 */
typedef struct mgn_moeadrbf_data mgnp_moeadrbf_data;
typedef struct mgnp_moeadrbf_mdl_p mgnp_moeadrbf_mdl_p; // moead model params

// TODO add function pointers as types
typedef void (*mgn_kernel_f)(gsl_vector *r, double s);


struct mgnp_features_moeadrbf {
    gsl_vector *sigma;
};

// pointer passing helper
typedef struct mgn_pointer {
    void* p;
} mgn_ptr;

struct mgnp_rbf_weigts {
    gsl_vector *p_s;
    gsl_matrix *p_m_w; // used for internal model
    gsl_matrix *m_phi;
};

struct mgn_moeadrbf_data {
    pmgn_moeadrbf_templatep()
    /* private */
    size_t max_run;
    mgn_pop_matrix *tset; // tset pop size(nt)
    mgnLimit *search_lim; // for training pop init
    bool l_lhci_lim; // true to free lhci limits
//    gsl_matrix *pm;    // pop for eval with model size(n)
    int scalarf;
    size_t mdl_size;// scalarization function, this is fixed to PBI
    size_t mdl_k;
    size_t mdl_wsize;
    mgn_model_rbfdata *rbf_data;
    mgn_pop *p_aprox;
    mgn_count_ciclic sel_idx;
    mgnp_moeadrbf_mdl_p *model_data;
    mgn_popl *arc;
};

#include "mgn_mop.h"
#include "mgn_de.h"
#include "mgn_cluster_m.h"

// moead de min optim
struct mgnp_moeadrbf_s_optim_ef_p {
    mgnp_de_ef_param_m()
    gsl_matrix *y_t; //training matrix
//    gsl_vector *y_cur;
};

struct mgnp_moeadrbf_s_optim_p {
    mgn_mop_param_common();
    mgn_pop_matrix *tset; //training matrix
    cluster_data *km;
    gsl_matrix *f_p;
    gsl_matrix *m_w;
    gsl_matrix *m_phi;
    mgn_kernel_f rbf;
};

// TODO add all dynamic parameters
// ======== MOEA/D model =========
struct mgnp_moeadrbf_mdl_p {
    size_t w_size;
    size_t mdl_size;
    size_t mdl_k; // scalarization function, this is fixed to PBI
    gsl_matrix *iW;
    mgnLimit *search_lim;
    mgn_model_rbfdata *rbf_data;
    gsl_vector *lambda;
    cluster_data *km;
    gsl_matrix *mphi;
    mgnLimit *mop_limits;
    size_t neighbours;
    size_t runs;
    mgn_ga_sets *ga_probs;
};



/*
 * Public Functions
 */

//common functions
mgn_pop_proto* mgn_moeadrbf_pop_get(mgnMoa *moa);

mgnp_moeadrbf_data* mgn_moeadrbf_features(mgnMoa* moa);

bool mgnp_moeadrbf_stop(mgnMoa *moa);

void mgnp_moeadrbf_set_ga_vals(mgnMoa *moa, mgn_ga_sets *ga);


/*
 * Private Functions
 */

/*
 * RBF optim S helpers
 */
int mgnp_moeadrbf_s_optim_mop_min(gsl_vector *f1
                                  ,gsl_vector *f2
                                  ,mgn_de_ef_param *ef_prm);

// x is sigma
// f is trained values y_predicted
// g is none
void mgnp_moeadrbf_s_optim_mop(gsl_vector *x
                               ,gsl_vector *f
                               ,gsl_vector *g
                               ,struct mgnp_moeadrbf_s_optim_p *prm
);

void mgnp_moeadrbf_optim_s(
    struct mgnp_rbf_weigts *vals
    , cluster_data *km
    , mgn_pop_matrix *tset
    , mgn_kernel_f rbf
);

// find the vector weidths lambda
// this could be done solving Ax=b system
gsl_vector *mgnp_moeadrbf_find_lambda(mgn_pop_matrix *tset
                                      , mgn_model_rbfdata *res
                                      , size_t size);

/*
 *
 * MOEAD model start
 *
 */
void pmgn_moeadrbf_lhci_update(mgnp_moeadrbf_data *rbfdata, mgn_pop_matrix *popm);

void mgnp_moeadrbf_mdl_defparam(mgnp_moeadrbf_data *moeadrbf
                                , size_t w_size
                                , double cr
                                , double mr);


void mgnp_moeadrbf_mdl_mop(gsl_vector *x
                           , gsl_vector *f
                           , gsl_vector *g
                           , mgnp_moeadrbf_mdl_p* prm);


void mgnp_moeadrbf_mdl_mop_param_free(mgnp_moeadrbf_data *moeadrbf);

void mgnp_moeadrbf_mdl_optim(mgnp_moeadrbf_mdl_p *mdl_p, mgn_pop *p_e, mgn_ptr *m_wl);


/*
 *
 * Selection and update
 *
 */
#include "uthash.h"

typedef struct mgnp_moeadrbf_selhash {
    int id;
    double* key;
    UT_hash_handle hh;
} mgnp_selhash;

// m_w matrix of w vectors in moead
// m_ws mgn_ptr with matrix of w vectors of moeadRBF
// sel: mgn_pop selected pop size = m_ws row size.
void mgnp_moeadrbf_select(mgn_count_ciclic *sel_idx
                          ,gsl_matrix *m_w
                          ,gsl_matrix *m_ws
                          ,mgn_pop *p_dest
                          ,mgn_pop* p_orig
                          );

void mgnp_moeadrbf_update_refine(mgn_pop *pop_newt, mgn_pop *pop_sel);

void mgnp_moeadrbf_update(mgnp_moeadrbf_data *mrbf, mgn_pop *pop_sel);

void mgnp_moeadrbf_pop_update(mgnp_moeadrbf_data *mrbf);



/*
 *
 * run common routines
 *
 */

struct mgnp_model_rbfdata {
    mgn_kernel_f kernel;
    struct mgnp_rbf_weigts mdl_rbf; // sigma, phi, w
};

// calculates phi and w values for rbf
// optimizes sigma
void pmgn_moeadrbf_calcNN(mgn_model_rbfdata *models
                          , mgn_pop_matrix *tset
                          , cluster_data *cdat
                          , cluster_data_extra *cdata_extra
                          );

void pmgn_moeadrbf_calcPhi(mgn_model_rbfdata *models
                           , mgn_pop_matrix *tset
                           , cluster_data *cdat
                           );


#endif //MOGEN_MGN_MOEAD_RBF_COMMON_H


