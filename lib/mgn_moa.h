#ifndef _LIB_MGN_MOA_H_
#define _LIB_MGN_MOA_H_

#include <stdbool.h>
#include <stdio.h>

#include "mgn_types.h"

#define MOA_NAME_LEN 32

typedef struct mgn_moa_t mgnMoa;
typedef struct mgn_moa_ga_set mgn_ga_sets;

struct mgn_moa_t {
    char name[MOA_NAME_LEN];
    size_t tot_exec;
    void (*run)(mgnMoa*);
    bool (*stop)();
    void (*set_ga_vals)(mgnMoa*, mgn_ga_sets*);
    void* features;
};

struct mgn_moa_ga_set {
    double cross_rate;
    double mut_rate;
    double *mut_llim;
    double *mut_ulim;
};

bool mgn_moa_solve(mgnMoa *moa, size_t runs)
{
    size_t i = 0;
    for (i = 0; i < runs; ++i) {
        moa->run(moa);
    }
    printf("total exec: %zu\n", moa->tot_exec);
    printf("gens: %zu\n", i);
    return true;
}

#endif // _LIB_MGN_MOA_H_
