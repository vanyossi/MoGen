
#include <stdio.h>

#include "mgn_random.h"

#include "population.h"
#include "individual.h"

void test_init_value(void* in, void* param);

int main()
{
    printf("test ran!\n");
    // rnd_initialize();
    for (size_t i = 0; i < 12; i++)
    {
        printf("rand %g\n", rnd_getUniform());
    }
    
    IndvParam param = {12,2,0};
    IndvOps *iops = mgn_IndvOps_init();
    

    MgnPop *pop = mgn_pop_alloc(100, (void*)iops, &param);
    mgnLimit ilim = {12.5,36.3};
    mgn_pop_init(pop, mgn_ind_init_rand, &ilim);

    // printf("indv 0 address %p\n", pop->I);
    // printf("indv 0 address %p\n", &pop->I[0]);
    // printf("indv 0 address %p\n", &(((Individual*)(pop->I))[0]) );
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

    MgnPop *pop1 = mgn_pop_alloc(5, (void*)iops,&param);
    MgnPop *pop2 = mgn_pop_alloc(15,(void*)iops,&param);

    double val = 2.34;
    mgn_pop_init(pop1,test_init_value, &val);
    val = 5;
    mgn_pop_init(pop2,test_init_value, &val);

    MgnPop *join = mgn_pop_join(pop1,pop2);

    printf("%.5f ", gsl_vector_get(mgn_indv_getx_vec(join,2),0));
    printf("%.5f ", gsl_vector_get(mgn_indv_getx_vec(join,10),0));
    printf("\n");

    mgn_pop_free(pop1);
    mgn_pop_free(pop2);
    mgn_pop_free(join);

    mgn_IndvOps_free(iops);
    return 0;
}

void test_init_value(void* in, void* param)
{
    double value = *(double*)param;
    Individual *ind = (Individual*) in;
    gsl_vector *x = ind->x;

    for (size_t i = 0; i < x->size; i++) {
        gsl_vector_set(x,i, value);
    }
}

