#ifndef _LIB_MGN_MOA_H_
#define _LIB_MGN_MOA_H_

#include "mgn_types.h"

#include <stdbool.h>

#define MOA_NAME_LEN 32

//typedef struct mgn_moa_t mgnMoa;
//typedef struct mgn_moa_ga_set mgn_ga_sets;

struct mgn_moa_t {
    char name[MOA_NAME_LEN];
    size_t tot_exec;
    void (*run)(mgnMoa*);
    bool (*stop)();
    void (*set_ga_vals)(mgnMoa*, mgn_ga_sets*);
    void* features;
    mgn_pop* (*pop_get)(mgnMoa*);
};

struct mgn_moa_ga_set {
    double cross_rate;
    double mut_rate;
    double *mut_llim;
    double *mut_ulim;
};

bool mgn_moa_solve(mgnMoa *moa, size_t runs);

void mgn_moa_set_mop();

#endif // _LIB_MGN_MOA_H_
