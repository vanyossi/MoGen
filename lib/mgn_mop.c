//
// Created by Iván Yossi on 29/07/22.
//
#include "mgn_mop.h"
#include "population.h"

mgnMop* mgn_mop_alloc()
{
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

size_t mgn_mop_eval_pop_index(mgnMop *mop, mgn_pop *pop, void *params, size_t pos, size_t size)
{
    if (pos + size > pop->size) return 0;

    for (size_t i = pos; i < (pos + size); i++) {
        pop->ops->eval(mop, mgn_pop_get(pop,i),params);
    }
    return size;
}

// void eval_real(void *inx, void* inf, void *ing){
//     castDouble(inx,x);
//     castDouble(inf,f);
//     castDouble(ing,g);
// }
