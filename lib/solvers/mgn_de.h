/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */


#ifndef MOGEN_MGN_DE_H
#define MOGEN_MGN_DE_H

#include "mgn_moa.h"

typedef struct mgnp_de_ef_param mgn_de_ef_param;
#define mgn_de_ef(ef) int (*ef)(void*, void*, mgn_de_ef_param*)
#define mgn_cast_de_ef(fname) (int (*)(void*, void*, mgn_de_ef_param*))fname

#define mgnp_de_ef_param_m() \
    size_t pos;

struct mgnp_de_ef_param {
    mgnp_de_ef_param_m()
};

mgnMoa* mgn_moa_de_alloc(size_t Np
                         ,void* iops
                         ,void* iparams
                         ,double factor
                         ,double cr
);

void mgn_de_init(mgnMoa* moa
                 ,void (*apply)(void*, void*)
                 ,void* apply_params
);

void mgn_moa_de_free(mgnMoa *de);

void mgn_de_setmop(mgnMoa *de, mgnMop *mop, mgn_de_ef(ef), mgn_de_ef_param* ef_p);

void mgn_de_eval(mgnMoa *de);

#endif //MOGEN_MGN_DE_H
