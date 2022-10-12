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

            moa->callback(moa);
        }
    }
    return true;
}

void mgn_moa_callback_virtual(mgnMoa* moa)
{
    UNUSED(moa);
}

mgn_pop_proto* mgn_moa_popget_virtual(mgnMoa* moa)
{
    UNUSED(moa);
    return 0;
}

mgnMoa* mgn_moa_alloc()
{
    mgnMoa *moa = malloc(sizeof(*moa));
    moa->callback = mgn_moa_callback_virtual;
    moa->pop_get = mgn_moa_popget_virtual;
    return moa;
}

void mgn_moa_set_mop(mgnMoa *moa, mgnMop *mop)
{
    moa->mop = mop;
}

void mgn_moa_set_callback(mgnMoa *moa, mgn_moa_callback_f callback)
{
    moa->callback = callback;
}
