#ifndef _MGN_POPULATION_
#define _MGN_POPULATION_

#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include "mgn_types.h"

enum _mgn_pop_type {
    MGN_REAL =          1 << 0,     // 0b000000001
    MGN_INT =           1 << 1,     // 0b000000010
    MGN_BIN =           1 << 2,     // 0b000000100
    MGN_MIX =           8 << 0,   // 0b000001000
};

struct _mgn_i_ops {
    void (*alloc)(void*, void*, void*);
    void (*copy)(void*, void*);
    void (*free)(void*);
    size_t (*sizeofp)();
    void* (*get_iparams)(void*);
    void* (*get_iops)(void*);
    void (*eval)(mgnMop*, void*);
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


MgnPop* mgn_pop_alloc(size_t size, void*(*indv_ops)(void*), void *params)
{
    MgnPop* pop = calloc(1, sizeof(*pop));
    pop->ops = (struct _mgn_i_ops*)indv_ops;
    pop->size = size;
    pop->current = 0;
    pop->I = calloc(size, pop->ops->sizeofp());

    // printf("sizeof %zu\n", pop->ops->sizeofp());

    char *pin = (char*)pop->I;
    for (size_t i = 0; i < size; i++) {
        pop->ops->alloc((void*)pin,indv_ops,params);
        pin += pop->ops->sizeofp();
    }
    
    return pop;
}

void mgn_pop_free(MgnPop *pop)
{
    char *pin = (char*)pop->I;
    for (size_t i = 0; i < pop->size; i++) {
        pop->ops->free((void*)pin);
        // pop->ops->free(&pop->I[i]);
        pin += pop->ops->sizeofp();
    }
    free(pop->I);
    free(pop);
    return;
}

void* mgn_pop_get(MgnPop *pop, size_t pos)
{
    char *pin = (char*)pop->I;
    for (size_t i = 0; i < pos; i++) {
        pin += pop->ops->sizeofp();
    }
    return pin;
}

void mgn_pop_init(MgnPop *pop, void (*apply)(void*, void*), void* params)
{
    char *pin = (char*)pop->I;
    for (size_t i = 0; i < pop->size; i++) {
        apply((void*)pin, params);
        pin += pop->ops->sizeofp();
    }
}

// both population are initiated
void mgn_pop_copy(MgnPop *dest, MgnPop *orig, size_t destpos, size_t origpos, size_t size)
{
    if (dest->ops->sizeofp() != orig->ops->sizeofp()) { return; }

    size_t isize = dest->ops->sizeofp();
    char *pdest = (char*)dest->I + (isize * destpos);
    char *porig = (char*)orig->I + (isize * origpos);

    if ((destpos + size) > dest->size) { return;}

    for (size_t i = 0; (i < size) && (destpos + i < dest->size); ++i) {
        dest->ops->copy(pdest, porig);
        pdest += isize;
        porig += isize;
    }
}

MgnPop* mgn_pop_join(MgnPop *p1, MgnPop *p2)
{
    MgnPop* npop = mgn_pop_alloc(p1->size + p2->size,
                                 p1->ops->get_iops(p1->I),p1->ops->get_iparams(p1->I));
    mgn_pop_copy(npop,p1,0,0,p1->size);
    mgn_pop_copy(npop,p2,p1->size,0,p2->size);

    return npop;
}


#endif // _MGN_POPULATION_
