//
// Created by Iván Yossi on 24/04/22.
//

#include <stdbool.h>
#include <gsl/gsl_vector.h>

#include "mgn_types.h"

#ifndef MOGEN_MGN_POP_PROTO_H
#define MOGEN_MGN_POP_PROTO_H


#define cast_get_iparams(fname) (mgn_pop_param (*)(void*))fname
#define pop_get_iparam(POP,POS) POP->ops->get_iparams(mgn_pop_get(POP, POS))
#define pop_alloc_pop(SIZE, RPOP) mgn_pop_alloc(SIZE \
                            ,RPOP->ops->get_iops(RPOP->I)\
                            , RPOP->ops->get_iparams_pointer(RPOP->I))\

#define mgnt_pop_free(name) void (*name)(void*)
#define mgnt_pop_copy(name) void (*name)(void*, void*)
#define mgnt_pop_alloc(name) void* (*name)(void*, void*, void*)

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
    void* (*alloc)(void*, void*, void*);\
    void (*copy)(void*, void*);\
    void (*free)(void*);\
    size_t (*sizeofp)();               \
    void* (*get_iparams_pointer)(void*);\
    void* (*get_iops)(void*);\
    void (*eval)(mgnMop*, void*, void*); \
    mgn_pop_param (*get_iparams)(void*);\
    void* (*next)(void*);               \
    void* (*prev)(void*);                   \
    void (*set_next)(void*,void*);      \
    void (*set_prev)(void*,void*);      \
    int (*rank_sort)(const void*, const void*); \
    void (*set_iparams)(void*, mgn_pop_param);

struct _mgn_i_ops {
    mgn_pop_ops()
};
// TODO evaluate (goes in mop)
// copy pops

#endif //MOGEN_MGN_POP_PROTO_H
