/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#ifndef MOGEN_MGF_MOA_H
#define MOGEN_MGF_MOA_H

#include "mgf_global_types.h"
#include "mop_report.h"

struct moa_t;

struct moa_cm_container {
    double prob;
    double eta;
};

enum moa_stop_criterion_t {
    MGN_STOPIF_GEN,
    MGN_STOPIF_EXEC,
    MGN_STOP_SIZE_MARKER            // not used, only for extending enum
};

struct moa_type_t {
    int data_size;
    // interface
    void (*typealloc)(struct moa_t*);
    mbool (*stop)(struct moa_t* moa, MoaStopCriterion criterion);
    mbool (*run)(struct mop_t* mop);
    void (*free)(struct moa_t*);
    // traits: initialized null
    struct moa_cm_container (*get_crossover_vals)(Moa*);
    struct moa_cm_container (*get_mutation_vals)(Moa*);
};

struct moa_t {
    struct moa_type_t *type;
    char name[64];
    struct mop_t* mop;
    MopReportStats stops;
    char buffer_start;      //!< Custom Moa start address
};


struct moa_type_t *mgf_moatype_new(
    int data_size,
    void (*typealloc)(struct moa_t *),
    mbool (*stop)(struct moa_t* moa, MoaStopCriterion criterion),
    mbool (*run)(struct mop_t* mop),
    void (*free)(struct moa_t *));

struct moa_type_t* mgf_moa_std();

struct moa_t *mgf_moa_new(Mop *mop, char *name, struct moa_type_t *type);
//Moa *moa_init(struct mop_t *mop, char *name, MoaTypes type, size_t mem_size);

//void mgf_moa_init(struct moa_t *moa, Mop *mop);

struct moa_type_t* mgf_moa_type(struct moa_t* self);

mbool mgf_moa_run(Moa* moa);

void* mgf_moa_buffer(struct moa_t* self);

void mgf_moa_free_std(struct moa_t* moa);

void mgf_moa_free(struct moa_t* moa);



mbool moa_run(struct mop_t *mop);

void moa_stopat_gen(Moa *moa, unsigned int gen);

void moa_stopat_eval(Moa *moa, unsigned int eval);

mbool moa_stop(Moa *moa, MoaStopCriterion criterion);


#endif //MOGEN_MGF_MOA_H
