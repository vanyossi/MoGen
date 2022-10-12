#ifndef _LIB_MGN_MOA_H_
#define _LIB_MGN_MOA_H_

#include "mgn_types.h"

#include <stdbool.h>

#define MOA_NAME_LEN 32

//typedef struct mgn_moa_t mgnMoa;
//typedef struct mgn_moa_ga_set mgn_ga_sets;

typedef void (*mgn_moa_callback_f)(mgnMoa*);
typedef mgn_pop_proto* (*mgn_moa_popget_f)(mgnMoa*);

struct mgn_moa_t {
    char name[MOA_NAME_LEN];
    size_t c_run;
    size_t max_exec;
    size_t tot_exec;
    mgnMop *mop;
    void (*run)(mgnMoa*);
    bool (*stop)(mgnMoa*);
    void (*set_ga_vals)(mgnMoa*, mgn_ga_sets*);
    void* features;
    mgn_moa_popget_f pop_get;
    mgn_moa_callback_f callback;
};

struct mgn_moa_ga_set {
    double cross_rate;
    double mut_rate;
    double *mut_llim;
    double *mut_ulim;
    double pbm_n;
    double sbx_m;
};

mgnMoa* mgn_moa_alloc();

bool mgn_moa_solve(mgnMoa *moa, size_t runs);

void mgn_moa_set_mop(mgnMoa *moa, mgnMop *mop);

void mgn_moa_set_callback(mgnMoa *moa, mgn_moa_callback_f callback);

#endif // _LIB_MGN_MOA_H_
