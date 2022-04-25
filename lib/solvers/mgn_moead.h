#ifndef _MGN_MOEAD_SOLVER_
#define _MGN_MOEAD_SOLVER_

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <float.h>

#include "mgn_moa.h"
#include "population.h"
#include "mgn_poplist.h"
#include "individual.h"
#include "decomposition/mgn_weights.h"
#include "decomposition/mgn_scalarization.h"
#include "operators/mgn_vector_distance.h"
#include "operators/mgn_gen_operator.h"

// TODO integrate to crossover struct
#ifndef MGN_SBX_N
#define MGN_SBX_N 20
#endif

#ifndef MGN_PBM_N
#define MGN_PBM_N 20
#endif

typedef struct moead_features moeadf;

struct moead_features {
    bool isalloc;
    bool ismopset;
    bool isprobset;
    mgn_ga_sets *ga_set;
    size_t size_nei;
    double *z;
    mgn_pop *pop;
    mgn_pop *epop;
    mgnMop *mop;
    mgnMop *_mop; // private
    gsl_matrix *wei;
    gsl_matrix *dist;
    gsl_matrix_int *dindex;
    double (*scalarize)(gsl_vector*, gsl_vector*, gsl_vector*);
};

moeadf* mgn_moead_getfeatures(mgnMoa* moead);
mgn_pop* mgn_moead_getpop(mgnMoa* moead);
void mgnp_moead_update_z(gsl_vector* x, gsl_vector* f, gsl_vector* g, moeadf* param);



void moead_reproduction(moeadf* set, mgn_pop* y_pop, size_t Ni)
{
    // crossover
    gsl_vector_int_view cur = gsl_matrix_int_row(set->dindex,Ni);
//    gsl_vector_int_view cur = gsl_matrix_int_subrow(set->dindex,Ni,1,set->dindex->size1 -1);

    gsl_vector_int *rand_sel = gsl_vector_int_alloc(set->size_nei);

//    memcpy(rand_sel->data, cur.vector.data, sizeof(int) * set->size_nei);
//    gsl_ran_shuffle(rnd_get_generator(), rand_sel->data, rand_sel->size,sizeof(int));

    mgn_select_nrandom_el(&cur.vector,rand_sel);

    mgn_pop* l_pop = pop_alloc_pop(2,set->pop);
    //          ,set->pop->ops->get_iparams( mgn_pop_get(set->pop,gsl_vector_int_get(rand_sel,0))).x->data

    if(set->ga_set->cross_rate > rnd_getUniform()) {
        mgn_genop_sbx(MGN_SBX_N
           , pop_get_iparam(set->pop, rand_sel->data[0]).x->data
          ,set->pop->ops->get_iparams( mgn_pop_get(set->pop,gsl_vector_int_get(rand_sel,1))).x->data
          ,l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->data
          ,l_pop->ops->get_iparams(mgn_pop_get(l_pop,1)).x->data
          ,l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->size
          , set->ga_set);
//        printf("data(1,2) %g, ", l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->data[0]);
//        printf("%g\n", l_pop->ops->get_iparams(mgn_pop_get(l_pop,1)).x->data[0]);

    } else {
        mgn_pop_copy(l_pop,set->pop,0,gsl_vector_int_get(rand_sel,0),1);
        mgn_pop_copy(l_pop,set->pop,1,gsl_vector_int_get(rand_sel,1),1);
    }
    gsl_vector_int_free(rand_sel);

    // mut_rate
    mgn_genop_pbm(MGN_PBM_N, set->ga_set->mut_rate
        ,l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->data
        ,set->ga_set->mut_llim
        ,set->ga_set->mut_ulim
        ,set->pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->size
        );

//    mgn_pop* ret_pop = pop_alloc_pop(1,l_pop);
    mgn_pop_copy(y_pop,l_pop,Ni, rnd_getUniform_int(2),1);
    mgn_pop_free(l_pop);
    return;
}

// return selected child index
void moead_update_neighbour(mgn_pop *lpop, moeadf *set, size_t Ni)
{
    gsl_vector_int_view Bi = gsl_matrix_int_row(set->dindex,Ni);
    gsl_vector_view zview = gsl_vector_view_array(set->z,set->wei->size2);
//    gsl_vector *wperm = gsl_vector_alloc(set->wei->size2);

//    size_t csel;
    for (size_t i = 0; i < set->size_nei; ++i) {
        int id = gsl_vector_int_get(&Bi.vector,i);
//        csel = rnd_getUniform_int(2); // pick one child at random
        gsl_vector_view wcur = gsl_matrix_row(set->wei,id);

//        gsl_vector_memcpy(wperm, &wcur.vector);
//        gsl_ran_shuffle(rnd_get_generator(),wperm->data, wperm->size, sizeof(double));

        double g1 = set->scalarize(&wcur.vector
            , lpop->ops->get_iparams(mgn_pop_get(lpop,Ni)).f
            , &zview.vector
            );
        double g2 = set->scalarize(&wcur.vector
           , set->pop->ops->get_iparams(mgn_pop_get(set->pop,id)).f
           , &zview.vector
        );

        if (g1 < g2) {
            mgn_pop_copy(set->pop,lpop,id,Ni,1);
        }
    }
//    gsl_vector_free(wperm);
}

void moead_update_ep(moeadf *set, mgn_pop *lpop)
{
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

/*
    // TODO make this an option
//    int max_ep_size = 200;
    // workaround as pop is not linkedlist
    // create new pop with only nondominated
    // exchange pop->I pointers
    // avoids loosing epop original pointer
    void* ind = set->pop->I;
//    mgn_pop *tpop = mgn_pop_alloc(1,set->pop->ops->get_iops(ind), set->pop->ops->get_iparams_pointer(ind));
//    mgn_pop_copy(tpop,lpop,0,0,1);

    for (size_t i = 0; i < set->epop->size; ++i) {
        int loc_dom_glob = vector_dominate(pop_get_iparam(lpop,0).f, pop_get_iparam(set->epop,i).f);
        if (loc_dom_glob < 0) {
//            printf("idom %d\n", loc_dom_glob);
            mgn_pop_param param = {1,0,0,0,0};
            set->epop->ops->set_iparams(mgn_pop_get(set->epop,i),param);
        }
        if (pop_get_iparam(lpop,0).rank == 0) {
            int glob_dom_loc = vector_dominate(pop_get_iparam(set->epop,i).f,pop_get_iparam(lpop,0).f);
            if (glob_dom_loc < 0) {
//                printf("itheydom %d\n", glob_dom_loc);
                mgn_pop_param param = {1,0,0,0,0};
                lpop->ops->set_iparams(mgn_pop_get(lpop,0),param);
            }
        }
    }
    mgn_pop *jpop = mgn_pop_join(lpop,set->epop);
    mgn_pop_qsort(jpop, indv_rank_sort);
    // TODO no need to rank sort all pop
    // make one pass to remove (add rank +1)
    // in same pass check if I will be added
//    mgn_pop_prank_sort(jpop);
//    int imdom = 0;

    // count non dominated
    int nondom = 0;
    for (size_t i = 0; i < jpop->size; ++i) {
        int rank = jpop->ops->get_iparams(mgn_pop_get(jpop,i)).rank;
        if (rank > 0 || nondom > 2000) break; // max epop artificial limit
        nondom++;
    }
//    printf("nondom %d\n", nondom);
// TODO limit ep population size
    mgn_pop *non_dom_pop = mgn_pop_alloc(nondom,set->pop->ops->get_iops(ind), set->pop->ops->get_iparams_pointer(ind));
    mgn_pop_copy(non_dom_pop,jpop,0,0,nondom);
    mgn_pop_exchange_iarray(set->epop, non_dom_pop);

    mgn_pop_free(non_dom_pop);
    mgn_pop_free(jpop);
//    mgn_pop_free(tpop);
*/

}

//void moead_run(mgnMop *mop, void* features)
void moead_run(mgnMoa* moead)
{
    moeadf* feat = mgn_moead_getfeatures(moead);
    // need check since init alloc are external calls
    if (!(feat->isalloc || feat->ismopset || feat->isprobset)) {
        // set stop condition to true
        return;
    }

    size_t n_weights = feat->wei->size1;
    mgn_pop* ypop = pop_alloc_pop(n_weights, feat->pop);
    for (size_t i = 0; i < n_weights; ++i) {
        moead_reproduction(feat, ypop,i);
        // TODO add mop params to mop.
        // update z
        moead->tot_exec += mgn_mop_eval_pop_index(feat->mop,ypop,NULL,i,1);
        // update rank
        mgn_mop_eval_pop_index(feat->_mop, ypop, feat,i,1);

        // update neighbour
        moead_update_neighbour(ypop,feat,i);

    }
    moead_update_ep(feat, ypop);

    mgn_pop_free(ypop);

    /* general framework assuming Tcheby
     *  sea W un conjunto de N vectores de peso
     *  sea zmin el punto de referencia
     *
     *  moead minimiza en una sola corrida N problemas
     *
     *  se define un vecindario para W_i como un conjunto de
     *      de los vectore de peso más cercanos
     *
     *  Población las mejores soluciones hasta ahora.
     *  En cada gen tenemos
     *      |P| = N
     *      FV_i = Fitnes values of x^i
     *      z = el mejor valor encontrado para f_i
     *      PE = población externa para almacenar soluciones no dom
     *
     *
     *  --ALGO
     *  input(mop, stop, N subprob, W, T = num vectors en vecidario)
     *  1.1 EP = 0
     *  1.2 Calcula euclidean entre los vectores y encuentra los T mas
     *      cercanos.
     *  1.3 Genera la población inicial t obten FV
     *  1.4 init z
     *
     *  2 Actualiza
     *  2.1 Reproducción, select random j,k y reproducirlos -> y
     *  2.2 Mejora:, aplica heuristica de mejora para y -> y^i
     *  2.3 Actualiza z, por cada j en m, para y^i z_j = f(y^i)
     *  2.4 Actualiza el vecindario.
     *  2.5 Actualiza EP
     *      remover todos las solcuiones dom por f(y^i)
     *      alade f(y^i) si nadie le domina.
     */

    /*
     * calculamos la distancia de todos los vectores contra todos (matrix)
     * Para W_i obtenemos los T más cercanos.
     *
     */

}

bool moead_stop()
{
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
// xfg param
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

// rpop = external reference pop
moeadf* mgn_moead_alloc_features(size_t H, size_t nobj, size_t T, mgn_pop *rpop,mgnMop *mop)
{
    moeadf *fe = malloc(sizeof(*fe));
    fe->scalarize = mgn_scalar_tchebycheff;
    fe->isprobset = false;
    fe->ismopset = false;

    // generate weight, get dist and sort;
//    fe->wei = mgn_weight_slattice_perm(H,nobj);

    fe->wei = mgn_weight_slattice_comb(H,nobj);
    for (size_t i = 0; i < fe->wei->size1; ++i) {
        gsl_vector_view wcur = gsl_matrix_row(fe->wei,i);
        gsl_ran_shuffle(rnd_get_generator(),
                        wcur.vector.data, wcur.vector.size,
                        sizeof(double)
                        );
    }
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
    fe->epop = rpop;
    fe->z = mgnp_moead_alloc_z(nobj);
    fe->size_nei = T;

    fe->mop = mop;
    fe->_mop = mgn_mop_alloc();
    fe->_mop->eval = mgn_cast_eval(mgnp_moead_update_z);
    fe->isalloc = true;
    return fe;
}

void mgn_moead_free(mgnMoa *moead){
    moeadf* fe = mgn_moead_getfeatures(moead);
    gsl_matrix_free(fe->wei);
    gsl_matrix_free(fe->dist);
    gsl_matrix_int_free(fe->dindex);
    mgn_pop_free(fe->pop);
    free(fe->z);
    mgn_mop_free(fe->_mop);
    free(fe);
    free(moead);
}

void moead_set_prob(mgnMoa*moead, mgn_ga_sets *gasets){
    moeadf* feat = mgn_moead_getfeatures(moead);
    feat->ga_set = gasets;
    feat->isprobset = true;
}

void moead_pop_evaluate(mgn_pop *pop, moeadf* set)
{
    mgn_mop_eval_pop(set->mop, pop, set);
    mgn_mop_eval_pop(set->_mop, pop, set);
}

mgnMoa* mgn_moead_init(size_t H,size_t nobj, size_t T, mgn_pop *rpop, mgnMop *mop
    ,void (*apply)(void*, void*), void* params)
{
    mgnMoa* moead = (mgnMoa*)calloc(1, sizeof(mgnMoa));
    strncpy(moead->name, "MOEA/D", MOA_NAME_LEN);
    moead->tot_exec = 0;
    moead->run = moead_run;
    moead->stop = moead_stop;
    moead->set_ga_vals = moead_set_prob;
    moead->features = mgn_moead_alloc_features(H,nobj,T,rpop,mop);

    moeadf *fe = mgn_moead_getfeatures(moead);
    mgn_pop_init(fe->pop, apply,params); // internal pop
    //evaluate pop
    moead_pop_evaluate(rpop, fe); // external pop
    moead_pop_evaluate(fe->pop,fe);
    mgn_pop_prank_sort(fe->pop);

    return moead;
}

moeadf* mgn_moead_getfeatures(mgnMoa* moead)
{
    return (moeadf*)moead->features;
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

// initializes population and calculates initial z value
//void mgn_moead_pop_init(mgnMoa* moead, void (*apply)(void*, void*), void* params)
//{
//    moeadf *moead_sp = mgn_moead_getfeatures(moead);
//    mgn_pop *pop = moead_sp->pop;
//}


#endif // _MGN_MOEAD_SOLVER_
