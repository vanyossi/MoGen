#ifndef _MGN_INDIVIDUAL_
#define _MGN_INDIVIDUAL_

#include <stdbool.h>
#include <stdlib.h>
// #include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>

#include "mgn_random.h"
#include "mgn_types.h"
#include "mgn_mop.h"
#include "mgn_pop_init.h"
#include "mgn_pareto.h"

#define indvCallOp(self,operator) self->ops->operator(self)
#define indGetXSize(self) indvCallOp(self,getXSize)
#define indGetObjSize(self) indvCallOp(self,getObjSize)
#define indGetConsSize(self) indvCallOp(self,getConsSize)

typedef struct _mgn_indv_ops mgn_indv_ops;
typedef struct _mgn_indv mgn_indv;
typedef struct _mgn_indv_params mgn_indv_param;

struct _mgn_indv_params {
    size_t realSize;
    size_t objSize;
    size_t consSize;
};

struct _mgn_indv {
    int rank;
    bool feasable;
    mgn_indv_ops* ops;
    mgn_indv *next;
    mgn_indv *prev;
    struct _mgn_indv_params* params;
    /*unsigned int size[3];*/
    gsl_vector        *x;
    gsl_vector        *f;
    gsl_vector        *g;
};

// TODO make macro of mandatory ops
struct _mgn_indv_ops {
    //mandatory
    mgn_pop_ops()
    // interface for this Indtype
    size_t (*getXSize)(mgn_indv*);
    size_t (*getObjSize)(mgn_indv*);
    size_t (*getConsSize)(mgn_indv*);
};

// Helper operation functions
size_t getXSize(mgn_indv* ind) { return ind->params->realSize; }
size_t getObjSize(mgn_indv* ind) { return ind->params->objSize; }
size_t getConsSize(mgn_indv* ind) { return ind->params->consSize; }

// Mandatory functions
void * mgn_indv_alloc(void* indv, void* ops, void* params);
void mgn_indv_copy(void *into, void *infrom);
void mgn_indv_free(void *in);
size_t mgn_sizeofp();
void* mgn_indv_get_params_p(void*);
mgn_pop_param mgn_indv_get_params(mgn_indv *in);
void* mgn_indv_get_ops(void* indv);
void mgn_indv_eval(mgnMop *mop, void* indv, void* param);
void mgn_indv_setparams(void* indv, mgn_pop_param);
int indv_rank_sort(const void* indv_a, const void* indv_b);
void* mgn_indv_next(void *indv);
void mgn_indv_set_next(mgn_indv *dest, mgn_indv *in);
void* mgn_indv_prev(void *indv);
void mgn_indv_set_prev(mgn_indv *dest, mgn_indv *in);

mgn_indv_ops* mgn_indv_ops_init()
{
    mgn_indv_ops *iops = calloc(1,sizeof(*iops));
    iops->alloc = mgn_indv_alloc;
    iops->copy = mgn_indv_copy;
    iops->free = mgn_indv_free;
    iops->sizeofp = mgn_sizeofp;
    iops->get_iops = mgn_indv_get_ops;
    iops->get_iparams_pointer = mgn_indv_get_params_p;
    iops->eval = mgn_indv_eval;
    iops->get_iparams = cast_get_iparams(mgn_indv_get_params);
    iops->set_iparams = mgn_indv_setparams;
    iops->next = mgn_indv_next;
    iops->set_next = (void (*)(void*,void*))mgn_indv_set_next;
    iops->prev = mgn_indv_prev;
    iops->set_prev = (void (*)(void*,void*))mgn_indv_set_prev;
    iops->rank_sort = indv_rank_sort;

    return iops;
}

void* mgn_indv_next(void *indv)
{
    mgn_indv *ind = (mgn_indv*)indv;
    return ind->next;
}

void* mgn_indv_prev(void *indv)
{
    mgn_indv *ind = (mgn_indv*)indv;
    return ind->prev;
}

void mgn_indv_set_next(mgn_indv *dest, mgn_indv *in)
{
    if(dest != NULL) {
        dest->next = in;
    }
}

void mgn_indv_set_prev(mgn_indv *dest, mgn_indv *in)
{
    if(dest != NULL) {
        dest->prev = in;
    }
}

void mgn_indv_ops_free(mgn_indv_ops *iops)
{
    free(iops);
    return;
}

void * mgn_indv_alloc(void* indv, void* ops, void* params)
{
    if (indv == NULL) {
        indv = malloc(mgn_sizeofp());
    }
    mgn_indv *nindv = (mgn_indv*)indv;
    nindv->ops = (mgn_indv_ops*)ops;
    nindv->next = 0;
    nindv->prev = 0;

    struct _mgn_indv_params* param = (struct _mgn_indv_params*)params;

    nindv->rank = 0;
    nindv->params = param;
//    nindv->size[0] = param->realSize;
//    nindv->size[1] = param->objSize;
//    nindv->size[2] = param->consSize;

    nindv->x = gsl_vector_alloc(param->realSize);
    nindv->f = gsl_vector_alloc(param->objSize);
    nindv->g = gsl_vector_alloc(param->consSize);

    return nindv;
}

void mgn_indv_free(void *in)
{
    mgn_indv *nindv = (mgn_indv*)in;
    gsl_vector_free(nindv->x);
    gsl_vector_free(nindv->f);
    gsl_vector_free(nindv->g);

    return;
}

void mgn_indv_free_all(void* in)
{
    mgn_indv_free(in);
    free(in);
}

void mgn_indv_copy(void *into, void *infrom)
{
    mgn_indv *from = (mgn_indv*)infrom;
    mgn_indv *to = (mgn_indv*)into;

    to->ops = from->ops;
    to->feasable = from->feasable;
    to->rank = from->rank;

    to->params = from->params;

//    to->size[0] = indGetXSize(from);
//    to->size[1] = indGetObjSize(from);
//    to->size[2] = indGetConsSize(from);

    gsl_vector_memcpy(to->x, from->x);
    gsl_vector_memcpy(to->f, from->f);
    gsl_vector_memcpy(to->g, from->g);

//    cblas_dcopy((int)to->x->size, from->x->data,1, to->x->data,1);
//    cblas_dcopy((int)to->f->size, from->f->data,1, to->f->data,1);
//    cblas_dcopy((int)to->g->size, from->g->data,1, to->g->data,1);

//    gsl_blas_dcopy(from->x,to->x);
//    gsl_blas_dcopy(from->f,to->f);
//    gsl_blas_dcopy(from->g,to->g);

    return;
}

size_t mgn_sizeofp()
{
    return sizeof(mgn_indv);
}

void mgn_indv_eval(mgnMop *mop, void* indv, void* param)
{
    mgn_indv *in = (mgn_indv*)indv;
    if (mop->eval) {
        mop->eval(in->x, in->f, in->g, param);
    }
    return;
}

void* mgn_indv_get_ops(void* indv)
{
    mgn_indv *in = (mgn_indv*)indv;
    return in->ops;
}

void* mgn_indv_get_params_p(void* indv)
{
    mgn_indv *in = (mgn_indv*)indv;
    return in->params;
}

// indv helpers
gsl_vector* mgn_indv_getx_vec(mgn_pop *pop, size_t in)
{
    return ((mgn_indv *)mgn_pop_get(pop,in))->x;
}

gsl_vector* mgn_indv_geto_vec(mgn_pop *pop, size_t in)
{
    return ((mgn_indv *)mgn_pop_get(pop,in))->f;
}

gsl_vector* mgn_indv_getc_vec(mgn_pop *pop, size_t in)
{
    return ((mgn_indv *)mgn_pop_get(pop,in))->g;
}

mgn_indv* mgn_indv_get(mgn_pop *pop, size_t in)
{
    return (mgn_indv *)mgn_pop_get(pop,in);
}

void mgn_ind_init(void* in, void* none)
{
    UNUSED(none);
    mgn_indv *ind = (mgn_indv*) in;
    for (size_t i = 0; i < getXSize(ind); i++) {
        gsl_vector_set(ind->x, i, rnd_getUniform());
    }
    return;
}

void mgn_ind_init_rand(void* in, void* limits)
{
    mgn_indv *ind = (mgn_indv*) in;
    mgnLimit* lim = (mgnLimit*) limits;

    for (size_t i = 0; i < ind->x->size; i++) {
        gsl_vector_set(ind->x,i, rnd_getUniform_limit(lim->min,lim->max));
    }
}

// for qsort
int mgn_ind_pareto_sort(const void* indv_a, const void* indv_b)
{
    mgn_indv *ia = (mgn_indv*) indv_a;
    mgn_indv *ib = (mgn_indv*) indv_b;

    int dom = vector_dominate(ia->f, ib->f);
    return dom;
}

gsl_matrix *mgn_ind_matrix(mgn_pop *pop)
{
    gsl_matrix *M = gsl_matrix_alloc(pop->size, pop->ops->get_iparams(pop->I).f->size);
    char* ind = pop->I;
    for (size_t i = 0; i < pop->size; ++i) {
        gsl_matrix_set_row(M,i,pop->ops->get_iparams(ind).f);
        ind += pop->ops->sizeofp();
    }
    return M;
}

int indv_rank_sort(const void* indv_a, const void* indv_b)
{
    mgn_indv *ia = (mgn_indv*) indv_a;
    mgn_indv *ib = (mgn_indv*) indv_b;

    return (ia->rank < ib->rank)? -1 : 1;
}

void pmgn_indv_setrank(mgn_indv *ind, size_t rank)
{
    ind->rank = rank;
}

void mgn_pop_ind_setrank(void * indin, int rank)
{
    mgn_indv *ind = (mgn_indv*) indin;
    ind->rank = rank;
}

void mgn_indv_setparams(void* indin, mgn_pop_param param)
{
    mgn_indv *ind = (mgn_indv*) indin;
    if (param.rank >= 0) {
        pmgn_indv_setrank(ind,param.rank);
    }
}

void mgn_pop_prank_sort(mgn_pop *pop)
{
    gsl_matrix *M = mgn_ind_matrix(pop);
    int *dranks = gsl_matrix_pareto_rank(M); //alloc
//    mgn_indv *ind = pop->I;
    char* ind = pop->I;
    mgn_pop_param iparam;
    for (size_t i = 0; i < pop->size; ++i) {
        iparam.rank = dranks[i];
        pop->ops->set_iparams(ind,iparam);
//        mgn_pop_ind_setrank(ind,dranks[i]);
//        ind[i].rank = dranks[i];
        ind += pop->ops->sizeofp();
    }
    // TODO make prototype indv to add default actions
    mgn_pop_qsort(pop, pop->ops->rank_sort);
    free(dranks);
    gsl_matrix_free(M);
}

mgn_pop_param mgn_indv_get_params(mgn_indv *in)
{
    mgn_pop_param param;
    param.rank = in->rank;
    param.feasable = in->feasable;
    param.x = in->x;
    param.f = in->f;
    param.g = in->g;

    return param;
}

#endif // _MGN_INDIVIDUAL_
