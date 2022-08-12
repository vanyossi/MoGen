#ifndef _MGN_INDIVIDUAL_
#define _MGN_INDIVIDUAL_


#include "mgn_types.h"

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "mgn_pop_proto.h"


typedef struct pmgn_indv_ops mgn_indv_ops;
typedef struct pmgn_indv mgn_indv;
typedef struct pmgn_indv_params mgn_indv_param;


struct pmgn_indv_params {
    size_t realSize;
    size_t objSize;
    size_t consSize;
};

struct pmgn_indv {
    int rank;
    bool feasable;
    mgn_indv_ops* ops;
    mgn_indv *next;
    mgn_indv *prev;
    struct pmgn_indv_params* params;
    /*unsigned int size[3];*/
    gsl_vector        *x;
    gsl_vector        *f;
    gsl_vector        *g;
};

struct pmgn_indv_ops {
    //mandatory
    mgn_i_ops()
    // interface for this Indtype
    size_t (*getXSize)(mgn_indv*);
    size_t (*getObjSize)(mgn_indv*);
    size_t (*getConsSize)(mgn_indv*);
};

#define indvCallOp(self,operator) self->ops->operator(self)
#define indGetXSize(self) indvCallOp(self,getXSize)
#define indGetObjSize(self) indvCallOp(self,getObjSize)
#define indGetConsSize(self) indvCallOp(self,getConsSize)


mgn_indv_ops* mgn_indv_ops_init();

void mgn_ind_init(void* in, void* none);

void mgn_ind_init_rand(void* in, void* limits);

mgn_indv* mgn_indv_get(mgn_pop *pop, size_t in);

void* mgn_indv_get_ops(void* indv);
void mgn_indv_ops_free(mgn_indv_ops *iops);

void* mgn_indv_get_params_p(void* indv);

// get vector components
gsl_vector* mgn_indv_getx_vec(mgn_pop *pop, size_t in);
gsl_vector* mgn_indv_geto_vec(mgn_pop *pop, size_t in);
gsl_vector* mgn_indv_getc_vec(mgn_pop *pop, size_t in);

// Mandatory functions
void *mgn_indv_alloc(void* indv, void* ops, void* params);
void mgn_indv_copy(void *into, void *infrom);
void mgn_indv_free(void *in);
size_t mgn_sizeofp();
void* mgn_indv_get_params_p(void*);
mgn_pop_param mgn_indv_get_params(mgn_indv *in);
void* mgn_indv_get_ops(void* indv);
void mgn_indv_eval(mgnMop *mop, void* indv, void* param);
void mgn_indv_setparams(void* indv, mgn_pop_param param);
int indv_rank_sort(const void* indv_a, const void* indv_b);
void* mgn_indv_next(void *indv);
void mgn_indv_set_next(mgn_indv *dest, mgn_indv *in);
void* mgn_indv_prev(void *indv);
void mgn_indv_set_prev(mgn_indv *dest, mgn_indv *in);

void mgn_pop_prank_sort(mgn_pop *pop);

gsl_matrix *mgn_ind_matrix_f(mgn_pop *pop);
gsl_matrix *mgn_ind_matrix_x(mgn_pop *pop);

#endif // _MGN_INDIVIDUAL_
