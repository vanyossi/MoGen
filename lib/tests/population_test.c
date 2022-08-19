
#include <stdio.h>
#include <math.h>

#include "mgn_test.h"
#include "mgn_random.h"

#include "population.h"
#include "individual.h"

#include "mgn_poplist.h"
#include "mgn_pareto.h"

static mgn_indv_param param = {12, 2, 2};
mgn_indv_ops *iops;

bool test_initialization();
bool test_join();
bool test_rank_sort();
bool test_pop_exchange();
bool test_pop_copy();
bool test_popl_ops();

void init_value_f(void* in, void* iparam);

int main()
{
    iops = mgn_indv_ops_init();
    rnd_getUniform();

    mgn_test("pop initialization", test_initialization);
    mgn_test("pop join", test_join);
    mgn_test("pareto rank test", test_rank_sort);
    mgn_test("pop exchange", test_pop_exchange);
    mgn_test("pop copy", test_pop_copy);
    mgn_test("popl operations", test_popl_ops);

    printf("tests done! %d/%d passed\n", tests_pass, tests_run);


    mgn_indv_ops_free(iops);
    return 0;
}


void init_value_f(void* in, void* iparam)
{
    double value = *(double*)iparam;
    mgn_indv *ind = (mgn_indv*) in;
    gsl_vector *x = ind->x;

    for (size_t i = 0; i < x->size; i++) {
        gsl_vector_set(x,i, value);
    }
}

bool pop_check_x_val(mgn_pop *pop, double val){
    bool ok = true;
    for (size_t i = 0; i < pop->size; ++i) {
        gsl_vector *x = pop->ops->get_iparams(mgn_pop_get(pop,i)).x;
        for (int j = 0; j < x->size; ++j) {
            ok &= fabs(gsl_vector_get(x, j) - val) < 1e-4;
        }
    }
    return ok;
}

bool test_initialization()
{
    bool ok = true;
    mgn_pop *pop1 = mgn_pop_alloc(5, (void*)iops,&param);

    double val = 2.34;
    mgn_pop_init(pop1, init_value_f, &val);
    ok &= pop_check_x_val(pop1, val);

    val = 5;
    mgn_pop_init(pop1, init_value_f, &val);
    ok &= pop_check_x_val(pop1, val);

    mgn_pop_free(pop1);
    return ok;
}

bool test_join()
{
    bool ok = true;
    double val;

    mgn_pop *A = mgn_pop_alloc(5, (void*)iops,&param);
    val = 2;
    mgn_pop_init(A, init_value_f, &val);
    mgn_pop *B = mgn_pop_alloc(12, (void*)iops,&param);
    val = 3.4;
    mgn_pop_init(B, init_value_f, &val);

    mgn_pop *C = mgn_pop_join(A,B);

    ok &= C->size == 17;
    ok &= C->ops->get_iparams(mgn_pop_get(C,2)).x->data[0] == 2;
    ok &= C->ops->get_iparams(mgn_pop_get(C,15)).x->data[0] == 3.4;

    mgn_pop_free(A);
    mgn_pop_free(B);
    mgn_pop_free(C);
    return ok;
}

bool test_rank_sort()
{
    bool pass = true;
    int results[5] = {0, 0, 1, 2, 3};

    mgn_pop *pop = mgn_pop_alloc(5,(void*)iops,&param);

    double *f;
    // a=[0 1;2 5]; b=[3 9]; z=[1 0; 0 10]
    // a
    f = mgn_indv_geto_vec(pop,0)->data;
    f[0] = 2; f[1] = 5; // 2
    f = mgn_indv_geto_vec(pop,1)->data;
    f[0] = 3; f[1] = 9; // 3
    f = mgn_indv_geto_vec(pop,2)->data;
    f[0] = 1; f[1] = 0; // 0
    f = mgn_indv_geto_vec(pop,3)->data;
    f[0] = 0; f[1] = 10; // 1
    f = mgn_indv_geto_vec(pop,4)->data;
    f[0] = 0; f[1] = 1; // 0

    mgn_pop_prank_sort(pop);
    for (size_t i = 0; i < pop->size; ++i) {
        pass &= mgn_indv_get(pop,i)->rank == results[i];
    }

    mgn_pop_free(pop);
    return pass;
}

bool test_pop_exchange()
{
    bool pass = true;
    double val;

    mgn_pop *A = mgn_pop_alloc(5, (void*)iops,&param);
    mgn_pop *B = mgn_pop_alloc(8, (void*)iops,&param);

    val = 2; mgn_pop_init(A, init_value_f, &val);
    val = 8; mgn_pop_init(B, init_value_f, &val);

    mgn_pop_exchange_iarray(A,B);

    pass &= A->size == 8;
    pass &= B->size == 5;

    for (size_t i = 0; i < A->size; ++i) {
        pass &= pop_check_x_val(A, 8);
    }
    for (size_t i = 0; i < B->size; ++i) {
        pass &= pop_check_x_val(B, 2);
    }

    mgn_pop_free(A);
    mgn_pop_free(B);
    return pass;
}


bool test_pop_copy()
{
    bool pass = true;
    double val = 3;

    mgn_pop *A = mgn_pop_alloc(10, (void*)iops,&param);
    mgn_pop *B = mgn_pop_alloc(20, (void*)iops,&param);

    mgn_pop_init(A, mgn_ind_init, &param);
    mgn_pop_init(B, init_value_f, &val);

    for (int i = 0; i < A->size; ++i) {
        pass &= A->ops->get_iparams(mgn_pop_get(A, i)).x->data[0] != 3;
    }

    mgn_pop_copy(B,A,0,0,A->size);

    for (int i = 0; i < A->size; ++i) {
        pass &=
            B->ops->get_iparams(mgn_pop_get(B, i)).x->data[0]
            ==
            A->ops->get_iparams(mgn_pop_get(A, i)).x->data[0]
        ;
    }

    mgn_pop_free(A);
    mgn_pop_free(B);
    return pass;
}

bool test_popl_ops()
{
    // create popl
    // push 8 individuals
    bool pass = true;
    double val = 2.34;

    mgn_popl* popl = mgn_popl_alloc((void*)iops,&param);

    mgn_popl_push(
        popl, mgn_indv_alloc(NULL
                             , popl->ops->get_iops(popl->I)
                             , popl->ops->get_iparams_pointer(popl->I)));
    mgn_popl_push(popl, mgn_indv_alloc(NULL, (void *) iops, &param));

    for (int i = 0; i < 6; ++i) {
        mgn_popl_push(popl, mgn_indv_alloc(NULL, (void *) iops, &param));
    }

    pass &= popl->size == 8;

    // Set 4 to pos 1,2,6
    mgn_popl_cursor_start(popl);
    for (size_t i = 0; i < popl->size; ++i) {
        val = (double)i * 10;
        if (i == 1 || i == 2 || i == 6) {
            val = 4;
        }
        init_value_f(mgn_popl_next(popl), &val);
    }


    pass &= popl->ops->get_iparams(mgn_popl_get(popl,1)).x->data[0] == 4;
    pass &= popl->ops->get_iparams(mgn_popl_get(popl,2)).x->data[0] == 4;
    pass &= popl->ops->get_iparams(mgn_popl_get(popl,6)).x->data[0] == 4;

//    mgn_popl_cursor_reset(popl);
//    for (size_t i = 0; i < popl->size; ++i) {
//        double myval = popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0];
//        printf("val %zu, %.6f\n",i, myval);
//    }

    // remove values with 4
    mgn_popl_cursor_reset(popl);
    int j = 0;
    void* ind = mgn_popl_current(popl);
    while(ind != 0) {
        double myval = popl->ops->get_iparams(ind).x->data[0];
        if (myval == 4) {
//            printf("\teal %d, %.6f %d\n",j, myval, popl->size);
            ind = mgn_popl_pop_current(popl);
            double vals = popl->ops->get_iparams(ind).x->data[0];
//            printf("removed ind %p %g\n", ind, vals);
            mgn_indv_free(ind);
            free(ind);
            ind = mgn_popl_current(popl);

        } else {
//            printf("else\n");
            mgn_popl_next(popl);
            ind = mgn_popl_current(popl);
            j++;
        }
    }

    pass &= popl->size == 5;
    mgn_popl_cursor_reset(popl);
    for (int i = 0; i < 5; ++i) {
        pass &= popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0] != 4;
    }

    val = 34.5;
    mgn_popl_push(popl, mgn_indv_alloc(0,(void*)iops, &param));
    init_value_f(mgn_popl_get_last(popl), &val);

    mgn_popl_alloc_last(popl);
    popl->ops->copy(mgn_popl_get_last(popl), mgn_popl_get(popl,1));

    pass &= popl->size == 7;
//    mgn_popl_cursor_reset(popl);
//    printf("size %d\n", popl->size);
//    for (size_t i = 0; i < popl->size; ++i) {
//        double myval = popl->ops->get_iparams(mgn_popl_next(popl)).x->data[0];
//        printf("val %zu, %.6f\n",i, myval);
////        mgn_popl_next(popl);
//    }
    mgn_popl_cursor_reset(popl);
    mgn_indv *popped = mgn_popl_pop(popl,0);
    while(popped != 0) {
//        double myval = popl->ops->get_iparams(popped).x->data[0];
        popped = mgn_popl_pop(popl,0);
//        printf("vval, %.6f\n", myval);
    }

    pass &= popl->size == 0;
    mgn_popl_free(popl);

    return pass;
}

