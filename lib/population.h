#ifndef _MGN_POPULATION_
#define _MGN_POPULATION_

#include "mgn_pop_proto.h"

struct _mgn_pop {
    mgn_pop_ops()
    unsigned int current;
};

mgn_pop* mgn_pop_alloc(size_t size, void*(*indv_ops)(void*), void *params);
void mgn_pop_free(mgn_pop *pop);
void* mgn_pop_get(mgn_pop *pop, size_t pos);

void mgn_pop_init(mgn_pop *pop, void (*apply)(void*, void*), void* params);
void mgn_pop_copy(mgn_pop *dest, mgn_pop *orig, size_t destpos, size_t origpos, size_t size);
void mgn_pop_exchange_iarray(mgn_pop *pop_a, mgn_pop *pop_b);

mgn_pop* mgn_pop_join(mgn_pop *p1, mgn_pop *p2);
void mgn_pop_qsort(mgn_pop *pop, int (*sort)(const void*, const void*));
int mgnp_pop_sort_1d(const void* ia, const void* ib);
void pop_sort_1d(mgn_pop* pop);

#endif // _MGN_POPULATION_
