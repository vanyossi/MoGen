/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_moead_common.h"


moeadf* mgn_moead_getfeatures(mgnMoa* moead)
{
    return (moeadf*)moead->features;
}

void moead_update_ep(moeadf *set, mgn_pop *lpop)
{
    if(set->epop == NULL) {
        return;
    }
//    printf("update ep\n");
    void* in_sert = mgn_pop_get(lpop,0); // pop has size 1
    gsl_vector *in_fval = lpop->ops->get_iparams(in_sert).f;

    mgn_popl_cursor_reset(set->epop);
    while(mgn_popl_current(set->epop) != 0) {
        gsl_vector *ep_fval = set->epop->ops->get_iparams(mgn_popl_current(set->epop)).f;
        if(vector_dominate(in_fval,ep_fval) < 0) {
//            printf("dominated\n");
            void* indv = mgn_popl_pop_current(set->epop); // advances cursor
            set->epop->ops->free(indv);
            free(indv);
        } else {
//            printf("else");
            mgn_popl_next(set->epop);
        }
//        printf("size... %d .. .. %p\n", set->epop->size, mgn_popl_current(set->epop));
//        printf("updating...\n");
//        i++;
//        if (i > 10) {
//            break;
//        }
    }
    bool inset_in = true;
    mgn_popl_cursor_reset(set->epop);
//    printf("sizze pointer %p, %d", mgn_popl_current(set->epop), set->epop->size);
    while(mgn_popl_current(set->epop) != 0) {
//        printf("analyzing...\n");
        gsl_vector *ep_fval = set->epop->ops->get_iparams(mgn_popl_current(set->epop)).f;
        if(vector_dominate(ep_fval,in_fval) < 0) {
            inset_in = false;
            break;
        }
        mgn_popl_next(set->epop);
//        if (i > 10) { break; }
    }
    if (inset_in) {
//        printf("inser\n");
        void *last = mgn_popl_alloc_last(set->epop);
//        printf("pointer last %p\n", last);
        set->epop->ops->copy(last, in_sert);
//        printf("data last inserted %g\n", set->epop->ops->get_iparams(last).f->data[0]);
//        printf("ted!\n");
    }
    /*
    // for many solutions
    mgn_pop *jpop = mgn_pop_join(lpop,set->epop);
    mgn_pop_prank_sort(jpop);

    int nondom = 0;
    for (size_t i = 0; i < jpop->size; ++i) {
        int rank = jpop->ops->get_iparams(mgn_pop_get(jpop,i)).rank;
        if (rank > 0 || nondom > 10000) break; // max epop artificial limit
        nondom++;
    }

    mgn_pop *non_dom_pop = pop_alloc_pop(nondom,set->pop);
    mgn_pop_copy(non_dom_pop,jpop,0,0,nondom);
    mgn_pop_exchange_iarray(set->epop, non_dom_pop);

    mgn_pop_free(non_dom_pop);
    mgn_pop_free(jpop);
*/
}


bool moead_stop(mgnMoa *moa)
{
    UNUSED(moa);
    return false;
}


double* mgnp_moead_alloc_z(size_t size)
{
    double *z = calloc(size, sizeof(*z));
    for (size_t i = 0; i < size; ++i) {
        z[i] = DBL_MAX;
    }
    return z;
}


void mgnp_moead_update_z(gsl_vector* x, gsl_vector* f, gsl_vector* g, moeadf* param)
{
    UNUSED(x);
    UNUSED(g);
    double fi;

    for (size_t i = 0; i < f->size; ++i) {
        fi = gsl_vector_get(f,i);
        if (fi < param->z[i]) {
            param->z[i] = fi;
        }
    }
}


moeadf* mgn_moead_alloc_features(gsl_matrix *W
                                 , size_t nobj
                                 , size_t T
                                 , mgn_popl *rpop
                                 , mgnMop *mop
                                 , bool external)
{
    moeadf *fe = malloc(sizeof(*fe));
    fe->scalarize = mgn_scalar_tchebycheff;
    fe->isprobset = false;
    fe->ismopset = false;

    // generate weight, get dist and sort;
    fe->wei = W;

//    fe->wei = mgn_weight_slattice_comb(H,nobj);
//    for (size_t i = 0; i < fe->wei->size1; ++i) {
//        gsl_vector_view wcur = gsl_matrix_row(fe->wei,i);
//        gsl_ran_shuffle(rnd_get_generator(),
//                        wcur.vector.data, wcur.vector.size,
//                        sizeof(double)
//                        );
//    }

    fe->dist = gsl_vector_distance_matrix(fe->wei, 2.0);
    fe->dindex = gsl_matrix_int_alloc(fe->dist->size1, fe->dist->size2);

    for (size_t i = 0; i < fe->dist->size1; ++i) {
        gsl_vector_view crow = gsl_matrix_row(fe->dist,i);
        int *dorder = gsl_vector_qsort(&crow.vector);
        gsl_vector_int_view irank = gsl_vector_int_view_array(dorder,fe->dist->size2);
        gsl_matrix_int_set_row(fe->dindex,i,&irank.vector);
        free(dorder);
    }

    fe->pop = mgn_pop_alloc(fe->wei->size1,
        rpop->ops->get_iops(rpop->I),rpop->ops->get_iparams_pointer(rpop->I));
    fe->epop = (external)? rpop : NULL;
    fe->z = mgnp_moead_alloc_z(nobj);
    fe->size_nei = (fe->wei->size1 < T)? fe->wei->size1 : T;
//    printf("size wei %zu, %zu\n", fe->wei->size1, fe->size_nei);

    fe->mop = mop;
    fe->_mop = mgn_mop_alloc();
    fe->_mop->eval = mgn_cast_eval(mgnp_moead_update_z);
    fe->_mop->params = fe;
    fe->isalloc = true;
    return fe;
}


void mgn_moead_free(mgnMoa *moead)
{
    moeadf* fe = mgn_moead_getfeatures(moead);
//    gsl_matrix_free(fe->wei);
    gsl_matrix_free(fe->dist);
    gsl_matrix_int_free(fe->dindex);
    mgn_pop_free(fe->pop);
    free(fe->z);
    mgn_mop_free(fe->_mop);
    free(fe);
    free(moead);
}


void moead_set_prob(mgnMoa*moead, mgn_ga_sets *gasets)
{
    moeadf* feat = mgn_moead_getfeatures(moead);
    feat->ga_set = gasets;
    feat->isprobset = true;
}


void moead_pop_evaluate(mgn_pop *pop, moeadf* set)
{
    mgn_mop_eval_pop(set->mop, pop, set->mop->params);
    mgn_mop_eval_pop(set->_mop, pop, set);
}


mgnMoa* mgn_moead_common_init(gsl_matrix *W,size_t nobj, size_t T
                       ,mgn_popl *epop
                       ,mgnMop *mop
                       ,void (*apply)(void*, void*), void* params
                       ,bool external)
{
    mgnMoa* moead = mgn_moa_alloc();
    strncpy(moead->name, "MOEA/D", MOA_NAME_LEN);
    moead->tot_exec = 0;
//    moead->run = moead_run;
    moead->stop = moead_stop;
    moead->set_ga_vals = moead_set_prob;
    moead->features = mgn_moead_alloc_features(W,nobj,T,epop,mop,external);
    moead->mop = mop; // TODO pick one mop location

    moead->pop_get = mgn_moead_getpop;
//    mgn_pop_prank_sort(fe->pop);

    return moead;
}

void mgn_moead_pop_init_eval(mgnMoa *moa, mgn_initializer *init)
{
    moeadf *fe = mgn_moead_getfeatures(moa);
    if(fe->pop) {
        mgn_init_pop_lhc(fe->pop, init,0);

        moead_pop_evaluate(fe->pop,fe);
        moa->tot_exec += fe->pop->size;
    }
}

void mgn_moead_set_scalarization(mgnMoa* moead
                                 ,mgnf_decomp_scalar f)
{
    moeadf *moaf = mgn_moead_getfeatures(moead);
    moaf->scalarize = f;
}


gsl_matrix* mgn_moead_get_w(mgnMoa* moead)
{
    return ((moeadf*)moead->features)->wei;
}


mgn_pop* mgn_moead_getpop(mgnMoa* moead)
{
    return ((moeadf*)moead->features)->pop;
}

void mgn_moead_set_mop(mgnMoa* moead, mgnMop *mop, int type)
{
    UNUSED(type);
    moeadf *moaf = mgn_moead_getfeatures(moead);
    moaf->mop = mop;
}
