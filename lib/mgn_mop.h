// mgn_mop

/*
 * Mop does not need to know about pop
 * it only needs to know how to evaluate
 * and generate fitness
*/
#ifndef _MGN_MOP_
#define _MGN_MOP_

#include "mgn_types.h"

typedef struct mgnp_mop_param mgn_mop_param;

#define mgn_cast_eval(fname) (void (*)(void*, void*, void*, void*))fname
typedef void (*mgn_mop_f)(void*, void*, void*, void*);

#define castDouble(A,B) double* B = (double*)A

#define mgn_mop_param_common() \
    size_t pos;

typedef struct mgnp_mop_param {
    mgn_mop_param_common()
} mgn_mop_p;

struct _mgn_mop {
    char name[32];
    mgnLimit *limits;
    void* params;
    void (*eval)(void*, void*, void*,void*);
    void (*eval_array)(void*, void*, void*,void*);
    void (*free)(mgnMop*); // mandatory
};


mgnMop* mgn_mop_alloc();

void mgn_mop_free(mgnMop *mop);

size_t mgn_mop_eval_pop(mgnMop *mop, mgn_pop *pop, void *params);

size_t mgn_mop_eval_pop_index(mgnMop *mop, mgn_pop *pop,
                              void *params, size_t pos, size_t size);

#endif // _MGN_MOP_
