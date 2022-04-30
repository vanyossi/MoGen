
#include <stdio.h>

#include "mgn_random.h"

#include "population.h"
#include "individual.h"

#include "mgn_poplist.h"
//#include "mgn_pareto.h"

void test_init_value(void* in, void* param);

int main()
{
    printf("test ran!\n");
    // rnd_initialize();
    for (size_t i = 0; i < 12; i++)
    {
        printf("rand %g\n", rnd_getUniform());
    }
    
    mgn_indv_param param = {12, 2, 2};
    mgn_indv_ops *iops = mgn_indv_ops_init();


    mgn_pop *pop = mgn_pop_alloc(100, (void*)iops, &param);
    mgnLimit *ilim = mgn_limit_alloc(param.realSize);
    for (size_t i = 0; i < ilim->size; ++i) {
        ilim->min[i] = -1;
        ilim->max[i] = 1;
    }
    mgn_pop_init(pop, mgn_ind_init_rand, ilim);

    // printf("indv 0 address %p\n", pop->I);
    // printf("indv 0 address %p\n", &pop->I[0]);
    // printf("indv 0 address %p\n", &(((mgn_indv*)(pop->I))[0]) );
    // printf("aa %g\n", gsl_vector_get(
    //     mgn_indv_getx_vec_of(pop,10),
    //     0
    //     ));
    // gsl_vector_fprintf(stdout,mgn_indv_getx_vec(pop,10),"%.4f");
    for (size_t i = 0; i < pop->size; i++) {
        gsl_vector *tmp = mgn_indv_getx_vec(pop,i);
        printf("%.6f, %.6f\n",
        gsl_vector_min(tmp),gsl_vector_max(tmp));
    }
    mgn_pop_free(pop);

    mgn_pop *pop1 = mgn_pop_alloc(5, (void*)iops,&param);
    mgn_pop *pop2 = mgn_pop_alloc(15,(void*)iops,&param);

    double val = 2.34;
    mgn_pop_init(pop1,test_init_value, &val);
    val = 5;
    mgn_pop_init(pop2,test_init_value, &val);

    mgn_pop *join = mgn_pop_join(pop1,pop2);

    printf("%.5f ", gsl_vector_get(mgn_indv_getx_vec(join,2),0));
    printf("%.5f ", gsl_vector_get(mgn_indv_getx_vec(join,10),0));
    printf("\n");

    mgn_pop_free(pop1);
    mgn_pop_free(pop2);
    mgn_pop_free(join);

    double *f;
    pop = mgn_pop_alloc(4,(void*)iops,&param);
    // a=[0 1;2 5]; b=[3 9]; z=[1 0; 0 10]
    // a
    f = mgn_indv_geto_vec(pop,0)->data;
    f[0] = 2; f[1] = 5;
    f = mgn_indv_geto_vec(pop,1)->data;
    f[0] = 3; f[1] = 9;
    f = mgn_indv_geto_vec(pop,2)->data;
    f[0] = 1; f[1] = 0;
    f = mgn_indv_geto_vec(pop,3)->data;
    f[0] = 0; f[1] = 10;

//    mgnLimit nlim = {0,1};
//    mgn_pop_init(pop, mgn_ind_init_rand,&nlim);
    mgn_pop_prank_sort(pop);
    for (size_t i = 0; i < pop->size; ++i) {
        printf("d %d\n", mgn_indv_get(pop,i)->rank);
    }

//    gsl_matrix *i_matrix = mgn_ind_matrix(pop);
//    int *dom_index = gsl_matrix_pareto_rank(i_matrix);
//    for (size_t i = 0; i < pop->size; ++i) {
//        printf("d %d\n", dom_index[i]);
//    }
//    free(dom_index);
//    gsl_matrix_fprintf(stdout,i_matrix,"%.3f");
//    printf("Sizes %zu %zu", i_matrix->size1, i_matrix->size2);
//    gsl_matrix_free(i_matrix);

    mgn_pop_copy(pop,pop,0,3,1);
    printf("\n");
    printf("psize %d\n", pop->size);
    for (size_t i = 0; i < pop->size; ++i) {
        printf("in %zu: %d\n", i, mgn_indv_get(pop,i)->rank);
        gsl_vector_fprintf(stdout,mgn_indv_geto_vec(pop,i),"%.2f");
    }

    mgn_pop* exchange_pop = mgn_pop_alloc(0,(void*)iops,&param);;
    mgn_pop* newpop = mgn_pop_join(pop, exchange_pop);


//    printf("size %d %d", newpop->size, exchange_pop->size);
    for (size_t i = 0; i < newpop->size; ++i) {
        printf("aa %zu: %d\n", i, mgn_indv_get(newpop,i)->rank);
        gsl_vector_fprintf(stdout,mgn_indv_geto_vec(newpop,i),":: %.2f");
    }
    mgn_pop_exchange_iarray(pop,exchange_pop);

    mgn_pop_free(exchange_pop);
    mgn_pop_free(pop);
    mgn_pop_free(newpop);


    mgn_popl* popl = mgn_popl_alloc((void*)iops,&param);

//    mgn_popl_pop(popl,0);
    mgn_popl_push(
        popl, mgn_indv_alloc(NULL
                             , popl->ops->get_iops(popl->I)
                             , popl->ops->get_iparams_pointer(popl->I)));
    mgn_popl_push(popl, mgn_indv_alloc(0, (void *) iops, &param));

    for (int i = 0; i < 6; ++i) {
        mgn_popl_push(popl, mgn_indv_alloc(0, (void *) iops, &param));
    }
//    mgn_popl_push(popl, mgn_indv_alloc(0,(void*)iops, &param));
//    mgn_popl_push(popl, mgn_indv_alloc(0,(void*)iops, &param));
//
    mgn_popl_cursor_start(popl);
    printf("size %d\n", popl->size);
    for (size_t i = 0; i < popl->size; ++i) {
        val = i;
        if (i == 1 || i == 2 || i == 6) {
            val = 4;
        }
//        printf("pointer %p %p\n",popl->first, mgn_popl_next(popl));
        test_init_value(mgn_popl_next(popl),&val);
//        mgn_popl_next(popl);
    }
//    mgn_indv *tmp = mgn_popl_pop(popl,1);
//    mgn_indv_free(tmp);
//    free(tmp);
    mgn_popl_cursor_reset(popl);
    printf("size %d\n", popl->size);
    for (size_t i = 0; i < popl->size; ++i) {
        double myval = popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0];
        printf("val %zu, %.6f\n",i, myval);
//        mgn_popl_next(popl);
    }


    mgn_popl_cursor_reset(popl);
    int j = 0;
    void* ind = mgn_popl_current(popl);
    while(ind != 0) {
        double myval = popl->ops->get_iparams(ind).x->data[0];
        if (myval == 4) {
            printf("\teal %d, %.6f %d\n",j, myval, popl->size);
            ind = mgn_popl_pop_current(popl);
            double vals = popl->ops->get_iparams(ind).x->data[0];
            printf("removed ind %p %g\n", ind, vals);
            mgn_indv_free_all(ind);
            ind = mgn_popl_current(popl);

        } else {
            printf("else\n");
            mgn_popl_next(popl);
            ind = mgn_popl_current(popl);
            j++;
        }
    }
        mgn_popl_cursor_reset(popl);
    for (int i = 0; i < 5; ++i) {
        mgn_popl_next(popl);
    }
//    mgn_popl_current_pop(popl);

    mgn_popl_push(popl, mgn_indv_alloc(0,(void*)iops, &param));
//    val = 100;
//    test_init_value(mgn_popl_get(popl,12),&val);
    val = 34.5;
    test_init_value(mgn_popl_get_last(popl),&val);

    mgn_popl_cursor_reset(popl);
    printf("size %d\n", popl->size);
    for (size_t i = 0; i < popl->size; ++i) {
        double myval = popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0];
        printf("val %zu, %.6f\n",i, myval);
//        mgn_popl_next(popl);
    }
    mgn_popl_cursor_reset(popl);
    while(mgn_popl_current(popl) != 0) {
        double myval = popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0];
        printf("vval, %.6f\n", myval);
    }

    mgn_popl_free(popl);
    mgn_indv_ops_free(iops);
    return 0;
}

void test_init_value(void* in, void* param)
{
    double value = *(double*)param;
    mgn_indv *ind = (mgn_indv*) in;
    gsl_vector *x = ind->x;

    for (size_t i = 0; i < x->size; i++) {
        gsl_vector_set(x,i, value);
    }
}

