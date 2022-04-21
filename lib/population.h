#ifndef _MGN_POPULATION_
#define _MGN_POPULATION_

#include <stdlib.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "mgn_types.h"

#define cast_get_iparams(fname) (mgn_pop_param (*)(void*))fname
#define pop_get_iparam(POP,POS) POP->ops->get_iparams(mgn_pop_get(POP, POS))
#define pop_alloc_pop(SIZE, RPOP) mgn_pop_alloc(SIZE \
                            ,RPOP->ops->get_iops(RPOP->I)\
                            , RPOP->ops->get_iparams_pointer(RPOP->I))\

#define mgnt_pop_free(name) void (*name)(void*)
#define mgnt_pop_copy(name) void (*name)(void*, void*)
#define mgnt_pop_alloc(name) void (*name)(void*, void*, void*)
#define mgnt_pop_alloc(name) void (*name)(void*, void*, void*)

enum _mgn_pop_type {
    MGN_REAL =          1 << 0,     // 0b000000001
    MGN_INT =           1 << 1,     // 0b000000010
    MGN_BIN =           1 << 2,     // 0b000000100
    MGN_MIX =           8 << 0,   // 0b000001000
};

struct _mgn_pop_param_pointer {
    int rank;
    bool feasable;
    gsl_vector *x;
    gsl_vector *f;
    gsl_vector *g;
};

//TODO move ops to indv prototype

#define mgn_pop_ops() \
    void (*alloc)(void*, void*, void*);\
    void (*copy)(void*, void*);\
    void (*free)(void*);\
    size_t (*sizeofp)();               \
    void* (*get_iparams_pointer)(void*);\
    void* (*get_iops)(void*);\
    void (*eval)(mgnMop*, void*, void*); \
    mgn_pop_param (*get_iparams)(void*); \
    void (*set_iparams)(void*, mgn_pop_param);

struct _mgn_i_ops {
    mgn_pop_ops()
};
// TODO evaluate (goes in mop)
// copy pops

struct _mgn_pop {
    unsigned int size;
    unsigned int current;
    popType type;
    struct _mgn_i_ops *ops;
    void* I;
};


mgn_pop* mgn_pop_alloc(size_t size, void*(*indv_ops)(void*), void *params)
{
    mgn_pop* pop = calloc(1, sizeof(*pop));
    pop->ops = (struct _mgn_i_ops*)indv_ops;
    pop->size = size;
    pop->current = 0;
    pop->I = calloc(size, pop->ops->sizeofp());

    // printf("sizeof %zu\n", pop->ops->sizeofp());
    size_t isize = pop->ops->sizeofp();
    mgnt_pop_alloc(alloc) = pop->ops->alloc;
    char *pin = (char*)pop->I;
    for (size_t i = 0; i < size; i++) {
        alloc((void*)pin,indv_ops,params);
        pin += isize;
    }
    
    return pop;
}

void mgn_pop_free(mgn_pop *pop)
{
    char *pin = (char*)pop->I;
    size_t isize = pop->ops->sizeofp();
    mgnt_pop_free(ind_free) = pop->ops->free;
    for (size_t i = 0; i < pop->size; i++) {
        ind_free((void*)pin);
        // pop->ops->free(&pop->I[i]);
        pin += isize;
    }
    free(pop->I);
    free(pop);
    return;
}

void* mgn_pop_get(mgn_pop *pop, size_t pos)
{
    char *pin = (char*)pop->I;
    pin += pop->ops->sizeofp() * pos;
//    for (size_t i = 0; i < pos; i++) {
//        pin += pop->ops->sizeofp();
//    }
    return pin;
}

void mgn_pop_init(mgn_pop *pop, void (*apply)(void*, void*), void* params)
{
    char *pin = (char*)pop->I;
    size_t isize = pop->ops->sizeofp();
    for (size_t i = 0; i < pop->size; i++) {
        apply((void*)pin, params);
        pin += isize;
    }
}

// both population are initiated
void mgn_pop_copy(mgn_pop *dest, mgn_pop *orig, size_t destpos, size_t origpos, size_t size)
{
    if (dest->ops->sizeofp() != orig->ops->sizeofp()) { return; }

    size_t isize = dest->ops->sizeofp();
    char *pdest = (char*)dest->I + (isize * destpos);
    char *porig = (char*)orig->I + (isize * origpos);

    if ((destpos + size) > dest->size) { return;}

    mgnt_pop_copy(copy) = dest->ops->copy;
    for (size_t i = 0; (i < size) && (destpos + i < dest->size); ++i) {
        copy(pdest, porig);
        pdest += isize;
        porig += isize;
    }
}

void mgn_pop_exchange_iarray(mgn_pop *pop_a, mgn_pop *pop_b)
{
    void* I = pop_a->I;
    pop_a->I = pop_b->I;
    pop_b->I = I;

    int size;
    size = pop_a->size;
    pop_a->size = pop_b->size;
    pop_b->size = size;
}

mgn_pop* mgn_pop_join(mgn_pop *p1, mgn_pop *p2)
{
    mgn_pop* npop = mgn_pop_alloc(p1->size + p2->size,
                                 p1->ops->get_iops(p1->I),p1->ops->get_iparams_pointer(p1->I));
    mgn_pop_copy(npop,p1,0,0,p1->size);
    mgn_pop_copy(npop,p2,p1->size,0,p2->size);

    return npop;
}

void mgn_pop_qsort(mgn_pop *pop, int (*sort)(const void*, const void*))
{
    qsort(pop->I, pop->size, pop->ops->sizeofp(),sort);
}

#endif // _MGN_POPULATION_
