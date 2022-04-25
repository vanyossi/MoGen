//
// Created by Iván Yossi on 24/04/22.
//

#ifndef MOGEN_MGN_POPLIST_H
#define MOGEN_MGN_POPLIST_H

#include "mgn_types.h"
#include "mgn_pop_proto.h"

#include <stdlib.h>
#include <stdbool.h>

typedef struct mgn_popl_head mgn_popl;

struct mgn_popl_head {
    unsigned int size;
    popType type;
    struct _mgn_i_ops *ops;
    void *first;
    void *current;
    void* I; // for prototype creation
};

mgn_popl* mgn_popl_alloc(void*(*indv_ops)(void*), void *params)
{
    mgn_popl* pop = calloc(1, sizeof(*pop));
    pop->ops = (struct _mgn_i_ops*)indv_ops;
    pop->size = 0;
    pop->first = 0;
    pop->current = 0;

    mgnt_pop_alloc(alloc) = pop->ops->alloc;
    pop->I = alloc(0,indv_ops,params);

    return pop;
}

void* mgn_popl_get_last(mgn_popl *pop)
{
    void* current = pop->first;
    if(current != 0){
        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
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
        printf("lst %p %p\n", last, indv);
        pop->ops->set_next(last,indv);
        pop->ops->set_prev(indv,last);
    }
    pop->size++;

    return;
}

void mgn_popl_alloc_last(mgn_popl *pop){
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
}

void* mgn_popl_get(mgn_popl *pop, size_t pos)
{
    size_t ipos = 0;
    void* current = pop->first;
    if(current != 0){
        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
        while (ipos < pos  && current && ops->next(current) != 0) {
            current = ops->next(current);
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

void* mgn_popl_current(mgn_popl* popl){
    return popl->current;
}

void* mgn_popl_pop_current(mgn_popl* pop)
{
    void *popped = 0;
    if (pop->current != 0) {
        void* cu = pop->current;
        struct _mgn_i_ops *ops = pop->ops->get_iops(cu);
        popped = cu;
        ops->set_next(ops->prev(cu),ops->next(cu));
        ops->set_prev(ops->next(cu),ops->prev(cu));
        pop->current = ops->next(cu);
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
        struct _mgn_i_ops *ops = pop->ops->get_iops(pop->I);
        pop->current = ops->next(prev);
    }
    return prev;
}

void* mgn_popl_pop(mgn_popl *pop, size_t pos)
{
    size_t ipos = 0;
    void *popped = 0;
    if(pos < pop->size && pop->first != 0){
        void *current = pop->first;
        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
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
        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
        while (ops->next(current) != 0) {
            void* for_free = current;
            current = ops->next(current);
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

#endif //MOGEN_MGN_POPLIST_H
