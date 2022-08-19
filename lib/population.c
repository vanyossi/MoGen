//
// Created by Iván Yossi on 29/07/22.
//
#include "population.h"

#include <stdlib.h>

void* pmgn_pop_get(void* pop_in, size_t index);
void pmgn_pop_set(void* pop_in, void* indv, size_t index);

mgn_pop* mgn_pop_alloc(size_t size, void*(*indv_ops)(void*), void *params)
{
    mgn_pop* pop = calloc(1, sizeof(*pop));
    pop->ops = (struct mgn_i_ops*)indv_ops;
    pop->iparams = *((mgn_indv_param*)params);
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

    pop->get = pmgn_pop_get;
    pop->set = pmgn_pop_set;

    return pop;
}

void mgn_pop_free(mgn_pop *pop)
{
    if(pop->I != 0) {
        char *pin = (char*)pop->I;
        size_t isize = pop->ops->sizeofp();
        mgnt_pop_free(ind_free) = pop->ops->free;
        for (size_t i = 0; i < pop->size; i++) {
            ind_free((void*)pin);
            // pop->ops->free(&pop->I[i]);
            pin += isize;
        }
        free(pop->I);
    }
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

int mgnp_pop_sort_1d(const void* ia, const void* ib)
{
    struct mgn_iparam_container* a = (struct mgn_iparam_container*)ia;
    struct mgn_iparam_container* b = (struct mgn_iparam_container*)ib;

    return (a->val <= b->val)? -1 : 1;
}

void pop_sort_1d(mgn_pop* pop)
{
    struct mgn_iparam_container* fvals;
    fvals = calloc(pop->size, sizeof(*fvals));
    mgn_pop *tmp_pop = pop_alloc_pop(pop->size,pop);

    for (size_t i = 0; i < pop->size; ++i) {
        fvals[i].i = i;
        fvals[i].val = pop_get_iparam(pop,i).f->data[0];
    }
    qsort(fvals,pop->size, sizeof(*fvals), mgnp_pop_sort_1d);

    for (size_t i = 0; i < pop->size; ++i) {
        mgn_pop_copy(tmp_pop,pop,i,fvals[i].i,1);
    }

    mgn_pop_exchange_iarray(pop,tmp_pop);

    mgn_pop_free(tmp_pop);
    free(fvals);
}


void* pmgn_pop_get(void* pop_in, size_t index)
{
    mgn_pop *pop = (mgn_pop*)pop_in;
    return mgn_pop_get(pop,index);
}

void pmgn_pop_set(void* pop_in, void* indv, size_t index)
{
    mgn_pop *pop = (mgn_pop*)pop_in;
    pop->ops->copy(pop->get(pop,index),indv);
}

