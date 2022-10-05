//
// Created by Iván Yossi on 29/07/22.
//
#include "mgn_moa.h"

#include <stdio.h>


bool mgn_moa_solve(mgnMoa *moa, size_t runs)
{
    size_t i = 0;

    for (i = 0; i < runs; ++i) {
        if (moa->tot_exec <= moa->max_exec) {
            moa->run(moa);
            moa->c_run++;
        }
    }
    return true;
}

void mgn_moa_set_mop(mgnMoa *moa, mgnMop *mop)
{
    moa->mop = mop;
}
