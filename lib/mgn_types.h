#ifndef _MGN_TYPES_
#define _MGN_TYPES_

#include <stdlib.h>

#define checkFlag(data, flag) (data & flag) == flag
#define UNUSED(x) ((void)(x))

typedef struct pmgn_indv_params mgn_indv_param;

// individual parameters
struct pmgn_indv_params {
    size_t x_size;
    size_t f_size;
    size_t g_size;
};


typedef struct pmgn_pop_ops mgn_pop_proto;
typedef struct _mgn_pop mgn_pop;
typedef struct mgn_popl_head mgn_popl;

typedef struct _mgn_pop_param_pointer mgn_pop_param;
typedef enum emgn_pop_type popType;

typedef struct _mgn_mop mgnMop;
typedef struct mgn_moa_t mgnMoa;
typedef struct mgn_moa_ga_set mgn_ga_sets;

typedef struct mgn_limit mgnLimit;

struct mgn_limit {
    size_t size;
    double *min;
    double *max;
};

static mgnLimit* mgn_limit_alloc(size_t size)
{
    mgnLimit *limits = malloc(sizeof(*limits));
    limits->size = size;
    limits->min = calloc(size, sizeof(limits->min));
    limits->max = calloc(size, sizeof(limits->max));
    return limits;
}

static void mgn_limit_free(mgnLimit *limits)
{
    free(limits->min);
    free(limits->max);
    free(limits);
}



#endif // _MGN_TYPES_
