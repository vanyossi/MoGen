//
// Created by Iván Yossi on 29/07/22.
//




#include "mgn_moead.h"
#include "mgn_moead_common.h"


mgn_pop* moead_reproduction(moeadf* set, mgn_pop* y_pop, size_t Ni)
{
    UNUSED(y_pop);
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
        mgn_genop_sbx(set->ga_set->sbx_m
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
    mgn_genop_pbm(set->ga_set->pbm_n, set->ga_set->mut_rate
                  ,l_pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->data
                  ,set->ga_set->mut_llim
                  ,set->ga_set->mut_ulim
                  ,set->pop->ops->get_iparams(mgn_pop_get(l_pop,0)).x->size
    );

    mgn_pop* ret_pop = pop_alloc_pop(1,l_pop);
    mgn_pop_copy(ret_pop,l_pop,0, rnd_getUniform_int(2),1);
    mgn_pop_free(l_pop);
    return ret_pop;
}


// return selected child index
void moead_update_neighbour(mgn_pop *lpop, moeadf *set, size_t Ni)
{
//    printf("update nei\n");
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
            mgn_pop_copy(set->pop,lpop,id,0,1);
        }
    }
//    gsl_vector_free(wperm);
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
    mgn_pop* ypop = pop_alloc_pop(1, feat->pop);
    for (size_t i = 0; i < n_weights; ++i) {
        mgn_pop* l_pop = moead_reproduction(feat, ypop,i);
        // TODO add mop params to mop (?)
        moead->tot_exec += mgn_mop_eval_pop(feat->mop,l_pop,feat->mop->params);
        // update z
        mgn_mop_eval_pop(feat->_mop, l_pop, feat);
        // update rank

        // update neighbour
        moead_update_neighbour(l_pop,feat,i);

        moead_update_ep(feat, l_pop);
        mgn_pop_free(l_pop);
    }
    mgn_pop_free(ypop);

//    printf("freed ypop\n");

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

mgnMoa* mgn_moead_init(gsl_matrix *W,size_t nobj, size_t T
                       ,mgn_popl *epop
                       ,mgnMop *mop
                       ,void (*apply)(void*, void*), void* params
                       ,bool external)
{
    mgnMoa* moead = mgn_moead_common_init(W,nobj,T
                                          ,epop
                                          ,mop
                                          ,apply
                                          ,params
                                          ,external);
    moead->run = moead_run;

    return moead;
}















// initializes population and calculates initial z value
//void mgn_moead_pop_init(mgnMoa* moead, void (*apply)(void*, void*), void* params)
//{
//    moeadf *moead_sp = mgn_moead_getfeatures(moead);
//    mgn_pop *pop = moead_sp->pop;
//}
