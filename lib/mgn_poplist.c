//
// Created by Iván Yossi on 24/07/22.
//

#include "mgn_poplist.h"

#include <stdlib.h>
#include <stdbool.h>

void* pmgn_popl_get(void* pop_in, size_t index);
void pmgn_popl_set(void* pop_in, void* indv, size_t index);

mgn_popl* mgn_popl_alloc(void*(*indv_ops)(void*), void *params)
{
    mgn_popl* pop = calloc(1, sizeof(*pop));
    pop->ops = (struct mgn_i_ops*)indv_ops;
    pop->iparams = *((mgn_indv_param*)params);
    pop->size = 0;
    pop->first = 0;
    pop->current = 0;

    mgnt_pop_alloc(alloc) = pop->ops->alloc;
    pop->I = alloc(0,indv_ops,params);

    pop->get = pmgn_popl_get;
    pop->set = pmgn_popl_set;

    return pop;
}

void* mgn_popl_get_last(mgn_popl *pop)
{
    void* current = pop->first;
    if(current != 0){
        struct mgn_i_ops *ops = pop->ops->get_iops(current);
        while (current && ops->next(current) != 0) {
            current = ops->next(current);
        }
    }
    return current;
}

void mgn_popl_push(mgn_popl *pop, void *indv)
{
    void* last = mgn_popl_get_last(pop);
    if (last == 0) {
        pop->first = indv;
    } else {
//        printf("lst %p %p\n", last, indv);
        pop->ops->set_next(last,indv);
        pop->ops->set_prev(indv,last);
    }
    pop->size++;

    return;
}

// alloc last and returns pointer to last
void* mgn_popl_alloc_last(mgn_popl *pop)
{
    void* last = mgn_popl_get_last(pop);
//    printf("alloc last  %p\n", last);
    void* new = pop->ops->alloc(NULL
                                ,pop->ops->get_iops(pop->I)
                                ,pop->ops->get_iparams_pointer(pop->I));
//    printf("alloc new %p\n", new);

    if (last == 0) {
        pop->first = new;
    } else {
        pop->ops->set_next(last, new);
        pop->ops->set_prev(new, last);
    }
    pop->size++;
//    pop->last = new;

    return new;
}

void* mgn_popl_get(mgn_popl *pop, size_t pos)
{
    size_t ipos = 0;
    void* current = pop->first;
    if(current != 0){
//        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
        while (ipos < pos  && current && pop->ops->next(current) != 0) {
            current = pop->ops->next(current);
            ipos++;
        }
    }
    return current;
}

void mgn_popl_cursor_reset(mgn_popl *pop)
{
    pop->current = pop->first;
}

void mgn_popl_cursor_start(mgn_popl *pop)
{
    mgn_popl_cursor_reset(pop);
}

void* mgn_popl_current(mgn_popl* popl)
{
    return popl->current;
}

void* mgn_popl_pop_current(mgn_popl* pop)
{
    void *popped = 0;
    if (pop->current != 0) {
        void* cu = pop->current;
//        struct _mgn_i_ops *ops = pop->ops->get_iops(cu);
        popped = cu;
        pop->ops->set_next(pop->ops->prev(cu),pop->ops->next(cu));
        pop->ops->set_prev(pop->ops->next(cu),pop->ops->prev(cu));
        pop->current = pop->ops->next(cu);
//        printf("pop next %p\n", ops->next(cu));
        if (cu == pop->first) {
            pop->first = pop->current;
        }
    }
    pop->size--;
    return popped;
}

// return last current
void* mgn_popl_next(mgn_popl *pop)
{
    void* prev = pop->current;

    if(prev != 0) {
        //struct _mgn_i_ops *ops = pop->ops->get_iops(pop->I);
        pop->current = pop->ops->next(prev);
    }
    return prev;
}

void* mgn_popl_pop(mgn_popl *pop, size_t pos)
{
    struct mgn_i_ops *ops = pop->ops;

    size_t ipos = 0;
    void *popped = 0;
    if(pos < pop->size && pop->first != 0){
        void *current = pop->first;
//        struct mgn_i_ops *ops = pop->ops->get_iops(current);
        while (ipos < pos && current != 0) {
            current = ops->next(current);
            ipos++;
        }
        if(ipos == pos && current != 0) {
            popped = current;
            ops->set_next(ops->prev(current),ops->next(current));
            ops->set_prev(ops->next(current),ops->prev(current));

            if(pop->current == popped){
                pop->current = ops->next(popped);
            }
        }
        pop->size--;
    }
    return popped;
}

void mgn_popl_free(mgn_popl* pop)
{
    void* current = pop->first;
    mgnt_pop_free(ind_free) = pop->ops->free;
    if(current != 0){
//        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
        while (pop->ops->next(current) != 0) {
            void* for_free = current;
            current = pop->ops->next(current);
            ind_free(for_free);
            free(for_free);
        }
        ind_free(current);
        free(current);
    }
    ind_free(pop->I);
    free(pop->I);
    free(pop);

    return;
}

void* pmgn_popl_get(void* pop_in, size_t index)
{
    mgn_popl *popl = (mgn_popl*)pop_in;
    return mgn_popl_get(popl,index);
}

// TODO: should copy before setting
void pmgn_popl_set(void* pop_in, void* indv, size_t index)
{
    UNUSED(index);
    mgn_popl *popl = (mgn_popl*)pop_in;
    mgn_popl_push(popl,indv);
}
