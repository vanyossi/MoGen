
#include <stdio.h>

#include "mgn_random.h"

#include "population.h"
#include "individual.h"

#define UNUSED(x) ((void)(x))

void test_init_rand(void* in, void* param);

int main(int argc, char const *argv[])
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
    mgn_IndvOps_free(iops);

    return 0;
}

void test_init_rand(void* in, void* param)
{
    UNUSED(param);
    Individual *ind = (Individual*) in;
    gsl_vector *x = ind->x;

    double llim = 34;
    double ulim = 36;

    for (size_t i = 0; i < x->size; i++) {
        gsl_vector_set(x,i, rnd_getUniform_limit(llim,ulim));
    }

    return;
}

