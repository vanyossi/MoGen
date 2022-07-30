//
// Created by Iván Yossi on 24/04/22.
//

#ifndef MOGEN_MGN_POPLIST_H
#define MOGEN_MGN_POPLIST_H

#include "mgn_types.h"
#include "mgn_pop_proto.h"


//typedef struct mgn_popl_head mgn_popl;

struct mgn_popl_head {
    unsigned int size;
    popType type;
    struct _mgn_i_ops *ops;
    void *first;
    void *current;
    void* I; // for prototype creation
};

mgn_popl* mgn_popl_alloc(void*(*indv_ops)(void*), void *params);
void* mgn_popl_get_last(mgn_popl *pop);
void mgn_popl_push(mgn_popl *pop, void *indv);
void mgn_popl_alloc_last(mgn_popl *pop);
void* mgn_popl_get(mgn_popl *pop, size_t pos);

void mgn_popl_cursor_reset(mgn_popl *pop);
void mgn_popl_cursor_start(mgn_popl *pop);
void* mgn_popl_current(mgn_popl* popl);
void* mgn_popl_pop_current(mgn_popl* pop);
void* mgn_popl_next(mgn_popl *pop);
void* mgn_popl_pop(mgn_popl *pop, size_t pos);
void mgn_popl_free(mgn_popl* pop);

#endif //MOGEN_MGN_POPLIST_H
