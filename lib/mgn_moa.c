//
// Created by Iván Yossi on 29/07/22.
//
#include "mgn_moa.h"

#include <stdio.h>

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
