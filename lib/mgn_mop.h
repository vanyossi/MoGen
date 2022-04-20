// mgn_mop

/*
 * Mop does not need to know about pop
 * it only needs to know how to evaluate
 * and generate fitness
*/
#ifndef _MGN_MOP_
#define _MGN_MOP_

#include "mgn_types.h"
#include "population.h"

#define mgn_cast_eval(fname) (void (*)(void*, void*, void*, void*))fname

#define castDouble(A,B) double* B = (double*)A

struct _mgn_mop {
    char name[32];
    void (*eval)(void*, void*, void*,void*);
};

mgnMop* mgn_mop_alloc(){
    mgnMop* mop = calloc(1,sizeof(*mop));
    return mop;
}

void mgn_mop_free(mgnMop *mop)
{
    free(mop);
}

size_t mgn_mop_eval_pop(mgnMop *mop, mgn_pop *pop, void *params)
{
    for (size_t i = 0; i < pop->size; i++) {
        pop->ops->eval(mop, mgn_pop_get(pop,i),params);
    }
    return pop->size;
}

// void eval_real(void *inx, void* inf, void *ing){
//     castDouble(inx,x);
//     castDouble(inf,f);
//     castDouble(ing,g);
// }

#endif // _MGN_MOP_
