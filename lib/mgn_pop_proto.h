//
// Created by Iván Yossi on 24/04/22.
//

#ifndef MOGEN_MGN_POP_PROTO_H
#define MOGEN_MGN_POP_PROTO_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "mgn_types.h"


#define cast_get_iparams(fname) (mgn_pop_param (*)(void*))fname

#define pop_get_iparam(POP,POS) POP->ops->get_iparams(mgn_pop_get(POP, POS))
#define pop_get_iparamp(POP) POP->ops->get_iparams_pointer(POP->I)
#define pop_get_iops(POP) POP->ops->get_iops(POP->I)
#define pop_alloc_pop(SIZE, RPOP) mgn_pop_alloc(SIZE \
                            ,RPOP->ops->get_iops(RPOP->I)\
                            , RPOP->ops->get_iparams_pointer(RPOP->I))\

#define mgnt_pop_free(name) void (*name)(void*)
#define mgnt_pop_copy(name) void (*name)(void*, void*)
#define mgnt_pop_alloc(name) void* (*name)(void*, void*, void*)

enum emgn_pop_type {
    MGN_POP =          1 << 0,     // 0b000000001
    MGN_POPL =           1 << 1,     // 0b000000010
    MGN_BIN =           1 << 2,     // 0b000000100
    MGN_MIX =           8 << 0,   // 0b000001000
};

// TODO fix this naming mess
// individual parameters
struct pmgn_indv_params {
    size_t x_size;
    size_t f_size;
    size_t g_size;
};


struct _mgn_pop_param_pointer {
    int rank;
    bool feasable;
    gsl_vector *x;
    gsl_vector *f;
    gsl_vector *g;
};

// TODO rename all _structs to pstructs
// population operators
/*struct pmgn_pop_ops *met; \*/
#define mgn_pop_ops() \
    struct mgn_i_ops *ops; \
    unsigned int size; \
    void* I;          \
    mgn_indv_param iparams;        \
    void* (*get)(void*, size_t); \
    void (*set)(void*, void* indv, size_t idx);

struct pmgn_pop_ops {
    mgn_pop_ops()
};

// individual operators
//TODO move ops to indv prototype
#define mgn_i_ops() \
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

struct mgn_i_ops {
    mgn_i_ops()
};
// TODO evaluate (goes in mop)
// copy pops

struct mgn_iparam_container {
    size_t i;
    double val;
};


void mgn_pop_print(mgn_pop_proto *popp, FILE *stream);


#endif //MOGEN_MGN_POP_PROTO_H
