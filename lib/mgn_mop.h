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

#define castDouble(A,B) double* B = (double*)A

struct _mgn_mop {
    char name[32];
    void (*eval_vector_real)(gsl_vector*, gsl_vector*, gsl_vector*);
    void (*eval_real)(double*, double*, double*);
    void (*eval_integer)(void*, void*, void*);
    void (*eval_binary)(void*, void*, void*);
};

void mgn_mop_eval_pop(mgnMop *mop, MgnPop *pop)
{
    for (size_t i = 0; i < pop->size; i++) {
        pop->ops->eval(mop, mgn_pop_get(pop,i));
    }
}

// void eval_real(void *inx, void* inf, void *ing){
//     castDouble(inx,x);
//     castDouble(inf,f);
//     castDouble(ing,g);
// }

#endif // _MGN_MOP_