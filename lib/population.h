#ifndef _MGN_POPULATION_
#define _MGN_POPULATION_

#include <stdlib.h>
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
};

struct _mgn_pop {
    unsigned int size;
    unsigned int current;
    popType type;
    struct _mgn_i_ops *ops;
    void* I;
};


MgnPop* mgn_pop_alloc(int size, void*(*indv_ops)(void*), void *params)
{
    MgnPop* pop = calloc(1, sizeof(MgnPop));
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

#endif // _MGN_POPULATION_