#ifndef _MGN_TYPES_
#define _MGN_TYPES_

#include <stdlib.h>

#define checkFlag(data, flag) (data & flag) == flag
#define UNUSED(x) ((void)(x))

typedef struct _mgn_pop mgn_pop;
typedef struct _mgn_pop_param_pointer mgn_pop_param;
typedef enum _mgn_pop_type popType;

typedef struct _mgn_mop mgnMop;

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
