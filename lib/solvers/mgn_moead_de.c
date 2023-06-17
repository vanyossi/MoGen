/*
 *
 *  SPDX-FileCopyrightText: 2022 Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead_de.h"
#include "mgn_moead_common.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "mgn_gnuplot.h"

#define MAX_EPOP 2000

// DE operators
//same as de
void mgnp_moead_de_rndseq_except(gsl_vector_ulong *seq, size_t except)
{
    // generate seq
    size_t s_i = 0;
    for (size_t i = 0; i <= seq->size; ++i) {
        if(i == except) continue; // skip except element
        gsl_vector_ulong_set(seq,s_i,i);
        s_i++;
    }
    gsl_ran_shuffle(rnd_get_generator(), seq->data, seq->size, sizeof(seq->data[0]));
}


void mgn_adjust_boundaries_to(gsl_vector *dest, gsl_vector *ref, mgnLimit *limit)
{
    for (size_t i = 0; i < dest->size; ++i) {
        if (dest->data[i] < limit->min[i]) {
            dest->data[i] = limit->min[i]
                + rnd_getUniform() * (ref->data[i] - limit->min[i]);
            continue;
        }
        if (dest->data[i] > limit->max[i]) {
            dest->data[i] = limit->max[i]
                - rnd_getUniform() * (limit->max[i] - ref->data[i]);
        }
    }
}

mgn_pop* mgn_de_recombination(mgn_pop *pop
                          , size_t cur_i
                          , size_t T
                          , double cr
                          , double factor
                          , moeadf* set)
{
    // missing

    unsigned long jrand;

    gsl_vector_ulong *rindex = gsl_vector_ulong_alloc(pop->size - 1);
    mgnp_moead_de_rndseq_except(rindex,cur_i);
    gsl_vector_ulong_view rsel = gsl_vector_ulong_subvector(rindex,0,3);
    gsl_vector_ulong_set(&rsel.vector,2,cur_i);

    jrand = rnd_getUniform_int((int) pop->iparams.x_size);

    mgn_pop *c_pop = mgn_pop_alloc(1,pop->ops,&pop->iparams);
    // cross over needs: pop, rsel, jrand, params
    gsl_vector *r0 = pop->ops->get_iparams(mgn_pop_get(pop,rsel.vector.data[0])).x;
    gsl_vector *r1 = pop->ops->get_iparams(mgn_pop_get(pop,rsel.vector.data[1])).x;
    gsl_vector *r2 = pop->ops->get_iparams(mgn_pop_get(pop,rsel.vector.data[2])).x;
//    gsl_vector *xi = pop->ops->get_iparams(mgn_pop_get(pop,rsel.vector.data[3])).x;
    gsl_vector *tu = pop->ops->get_iparams(mgn_pop_get(c_pop,0)).x;

    for (size_t j = 0; j < tu->size; ++j) {
        if (rnd_getUniform() <= cr || j == jrand) {
            tu->data[j] = r2->data[j] + (factor * (r1->data[j] - r0->data[j]));
//            mgnp_de_mut(&fparam, j);
        } else {
            gsl_vector_set(tu,j,gsl_vector_get(r2,j));
        }
    }
    mgn_adjust_boundaries_to(tu,r2,set->mop->limits);

    gsl_vector_ulong_free(rindex);

    return c_pop;
}





// NOTE
// nr could be 0 if no de, and bigger than 0 for de
// for that we need a function to create the index arrays. from Bi
void moead_update_neighbour_de(mgn_pop *lpop
                               , moeadf *set
                               , size_t Ni
                               , size_t nr)
{

    gsl_vector_int_view Bi = gsl_matrix_int_row(set->dindex,Ni);
    gsl_vector_view zview = gsl_vector_view_array(set->z,set->wei->size2);

    gsl_permutation *perm = gsl_permutation_calloc(set->size_nei);
    gsl_ran_shuffle(rnd_get_generator(),perm->data,perm->size, sizeof(size_t));

    size_t csel = 0;
    for (size_t i = 0; i < set->size_nei && csel < nr; ++i) {

        int id = gsl_vector_int_get(&Bi.vector, gsl_permutation_get(perm,i));
//        csel = rnd_getUniform_int(2); // pick one child at random
        gsl_vector_view wcur = gsl_matrix_row(set->wei,id);

//        gsl_vector_memcpy(wperm, &wcur.vector);
//        gsl_ran_shuffle(rnd_get_generator(),wperm->data, wperm->size, sizeof(double));
        double theta = 0.5;
        double g1 = set->scalarize(&wcur.vector
                                   , lpop->ops->get_iparams(mgn_pop_get(lpop,0)).f
                                   , &zview.vector
                                   , &theta
        );
        double g2 = set->scalarize(&wcur.vector
                                   , set->pop->ops->get_iparams(mgn_pop_get(set->pop,id)).f
                                   , &zview.vector
                                   , &theta
        );

        if (g1 < g2) {
            csel++;
            mgn_pop_copy(set->pop,lpop,id,0,1);
        }
    }
    gsl_permutation_free(perm);
}




//void moead_run(mgnMop *mop, void* features)
void moead_de_run(mgnMoa* moead)
{
//    bool count = true; //tmp var while we count exec from moa
    moeadf* feat = mgn_moead_getfeatures(moead);
    // need check since init alloc are external calls
    if (!(feat->isalloc || feat->ismopset || feat->isprobset)) {
        // set stop condition to true
        return;
    }

    size_t n_weights = feat->wei->size1;
    mgn_pop* ypop = pop_alloc_pop(1, feat->pop);
    for (size_t i = 0; i < n_weights; ++i) {
//        if (moead->tot_exec >= moead->max_exec){
////            count = false;
//            break;
//        }
//        mgn_pop* l_pop = moead_de_reproduction(feat, ypop,i);
        mgn_pop* l_pop = mgn_de_recombination(feat->pop, i
                                              ,feat->size_nei
                                              ,feat->ga_set->cross_rate
                                              ,1.2
                                              ,feat);

        gsl_vector *indv_x = l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x;
        mgn_genop_pbm(feat->ga_set->pbm_n, feat->ga_set->mut_rate
                      ,indv_x->data
                      ,feat->ga_set->mut_llim
                      ,feat->ga_set->mut_ulim
                      ,indv_x->size
        );

        // TODO add mop params to mop (?)
        moead->tot_exec += mgn_mop_eval_pop(feat->mop,l_pop,feat->mop->params);
        // update z
        mgn_mop_eval_pop(feat->_mop, l_pop, feat);
        // update rank

        // update neighbour
        moead_update_neighbour_de(l_pop,feat,i,2);

        moead_update_ep(feat, l_pop);
        mgn_pop_free(l_pop);

        if(feat->epop && feat->epop->size > MAX_EPOP) {
            moead->tot_exec = moead->max_exec;
        }

//        if ( moead->tot_exec % 300 == 0) {
//            char *filename = malloc(64);
//            asprintf(&filename, "MOEAD_DE_%s_%zu", moead->mop->name, moead->tot_exec);
//            mgn_plot_fast((mgn_pop_proto*) feat->epop, filename, "sol");
//            free(filename);
//        }
    }

//    if (count) { moead->c_run++; }
    mgn_pop_free(ypop);
}

/**
 * Initializes moa parameters
 * @param epop external pop pointer
 * @param external, if TRUE training use external pop as training population
 *                  population initialization can be skipped
 */
mgnMoa* mgn_moead_de_init(gsl_matrix *W,size_t nobj, size_t T
                          ,mgn_popl *epop
                          ,mgnMop *mop
                          ,void (*apply)(void*, void*), void* params
                          ,bool external)
{
    mgnMoa* moead_de = mgn_moead_common_init(W,nobj,T
                        ,epop
                        ,mop
                        ,apply
                        ,params
                        ,external);
    strncpy(moead_de->name, "MOEAD-DE", MOA_NAME_LEN);
    moead_de->run = moead_de_run;
    moead_de->max_exec = 5000;
    moead_de->c_run = 0;
    moead_de->tot_exec = 0;

    return moead_de;
}
