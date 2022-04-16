#ifndef _MGN_INDIVIDUAL_
#define _MGN_INDIVIDUAL_

#include <stdbool.h>
#include <stdlib.h>
// #include <stdio.h>
#include <gsl/gsl_vector.h>

#include "mgn_random.h"
#include "mgn_types.h"
#include "mgn_mop.h"
#include "mgn_pop_init.h"


#define indvCallOp(self,operator) self->ops->operator(self)
#define indGetXSize(self) indvCallOp(self,getXSize)
#define indGetObjSize(self) indvCallOp(self,getObjSize)
#define indGetConsSize(self) indvCallOp(self,getConsSize)

typedef struct _mgn_indv_ops IndvOps;
typedef struct _mgn_indv Individual;
typedef struct _mgn_indv_params IndvParam;

struct _mgn_indv_params {
    int realSize;
    int objSize;
    int consSize;
};

struct _mgn_indv {
    int index;
    bool feasable;
    IndvOps* ops;
    struct _mgn_indv_params* params;
    /*unsigned int size[3];*/
    gsl_vector        *x;
    gsl_vector        *f;
    gsl_vector        *g;
};

struct _mgn_indv_ops {
    //mandatory
    void (*alloc)(void*, void*, void*);
    void (*copy)(void*, void*);
    void (*free)(void*);
    size_t (*sizeofp)();
    void* (*get_iparams)(void*);
    void* (*get_iops)(void*);
    void (*eval)(mgnMop*, void*);
    // interface for this Indtype
    size_t (*getXSize)(Individual*);
    size_t (*getObjSize)(Individual*);
    size_t (*getConsSize)(Individual*);
};

// Helper operation functions
int getRSize(Individual* ind) { return ind->params->realSize; }
int getObjSize(Individual* ind) { return ind->params->objSize; }
int getConsSize(Individual* ind) { return ind->params->consSize; }

// Mandatory functions
void mgn_indv_alloc(void* indv, void* ops, void* params);
void mgn_indv_copy(void *into, void *infrom);
void mgn_indv_free(void *in);
size_t mgn_sizeofp();
void* mgn_indv_get_params(void*);
void* mgn_indv_get_ops(void* indv);
void mgn_indv_eval(mgnMop *mop, void* indv);


IndvOps* mgn_IndvOps_init()
{
    IndvOps *iops = calloc(1,sizeof(*iops));
    iops->alloc = mgn_indv_alloc;
    iops->copy = mgn_indv_copy;
    iops->free = mgn_indv_free;
    iops->sizeofp = mgn_sizeofp;
    iops->get_iops = mgn_indv_get_ops;
    iops->get_iparams = mgn_indv_get_params;
    iops->eval = mgn_indv_eval;

    return iops;
}

void mgn_IndvOps_free(IndvOps *iops)
{
    free(iops);
    return;
}

void mgn_indv_alloc(void* indv, void* ops, void* params)
{
    Individual *nindv = (Individual*)indv;
    nindv->ops = (IndvOps*)ops;

    struct _mgn_indv_params* param = (struct _mgn_indv_params*)params;

    nindv->params = param;
//    nindv->size[0] = param->realSize;
//    nindv->size[1] = param->objSize;
//    nindv->size[2] = param->consSize;

    nindv->x = gsl_vector_alloc(param->realSize);
    nindv->f = gsl_vector_alloc(param->objSize);
    nindv->g = gsl_vector_alloc(param->consSize);

    return;
}

void mgn_indv_free(void *in)
{
    Individual *nindv = (Individual*)in;
    gsl_vector_free(nindv->x);
    gsl_vector_free(nindv->f);
    gsl_vector_free(nindv->g);

    return;
}

void mgn_indv_copy(void *into, void *infrom)
{
    Individual *from = (Individual*)infrom;
    Individual *to = (Individual*)into;

    to->ops = from->ops;
    to->feasable = from->feasable;

    to->params = from->params;

//    to->size[0] = indGetXSize(from);
//    to->size[1] = indGetObjSize(from);
//    to->size[2] = indGetConsSize(from);

    gsl_vector_memcpy(to->x, from->x);
    gsl_vector_memcpy(to->f, from->f);
    gsl_vector_memcpy(to->g, from->g);

    return;
}

size_t mgn_sizeofp()
{
    return sizeof(Individual);
}

void mgn_indv_eval(mgnMop *mop, void* indv)
{
    Individual *in = (Individual*)indv;
    if (mop->eval_vector_real) {
        mop->eval_vector_real(in->x, in->f, in->g);
    } else {
        mop->eval_real(in->x->data, in->f->data, in->g->data);
    }
    return;
}

void* mgn_indv_get_ops(void* indv)
{
    Individual *in = (Individual*)indv;
    return in->ops;
}

void* mgn_indv_get_params(void* indv)
{
    Individual *in = (Individual*)indv;
    return in->params;
}

// indv helpers
gsl_vector* mgn_indv_getx_vec(MgnPop *pop, size_t in)
{
    return ((Individual *)mgn_pop_get(pop,in))->x;
}

gsl_vector* mgn_indv_geto_vec(MgnPop *pop, size_t in)
{
    return ((Individual *)mgn_pop_get(pop,in))->f;
}

gsl_vector* mgn_indv_getc_vec(MgnPop *pop, size_t in)
{
    return ((Individual *)mgn_pop_get(pop,in))->g;
}

void mgn_ind_init(Individual* ind)
{
    for (size_t i = 0; i < indGetXSize(ind); i++) {
        gsl_vector_set(ind->x, i, rnd_getUniform());
    }
    return;
}

void mgn_ind_init_rand(void* in, void* limits)
{
    Individual *ind = (Individual*) in;
    mgnLimit* lim = (mgnLimit*) limits;

    for (size_t i = 0; i < ind->x->size; i++) {
        gsl_vector_set(ind->x,i, rnd_getUniform_limit(lim->min,lim->max));
    }
}


#endif // _MGN_INDIVIDUAL_
