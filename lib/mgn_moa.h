#ifndef _LIB_MGN_MOA_H_
#define _LIB_MGN_MOA_H_

#include <stdbool.h>

#include "mgn_types.h"

typedef struct _mgn_moa_t mgnMoa;

struct _mgn_moa_t {
    char name[32];
    void (*run)(mgnMop*);
    bool (*stop)();
    void* features; 
};

bool mgn_moa_solve(mgnMoa *moa, mgnMop *mop)
{
    moa->run(mop);
    return true;
}

#endif // _LIB_MGN_MOA_H_