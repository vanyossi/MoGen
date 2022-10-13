
#include <math.h>
#include <stdbool.h>

#include "indicadores.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_rstat.h"

#include "gsl_vector_additional.h"

static double m_p = 2;

// helpers
struct pind_func_dist {
    double (*f)(gsl_vector*, gsl_vector*, double);
    double pval;
};

// @TODO
// considerar solo soluciones no dominadas
// gsl_matrix_swap_rows <- para ordenar de no dom a dom.
// luego tomar sub_matrix
double pin_pow(double val, void *p)
{
    double inp = *(double*)p;
    return pow(val,inp);
}


gsl_matrix* repeatRowToSizeB(int row, gsl_matrix *A, gsl_matrix *B)
{
    gsl_matrix *Sc = gsl_matrix_alloc(B->size1, B->size2);
    gsl_matrix_set_all(Sc,1);
    gsl_vector *slice = gsl_vector_alloc(Sc->size2);
    for (size_t fp = 0; fp < Sc->size1; fp++) {
        gsl_vector_view cr = gsl_matrix_row(Sc,fp);
        gsl_matrix_get_row(slice,A,row);
        gsl_vector_mul(&cr.vector,slice);

    }
    gsl_vector_free(slice);
    return Sc;
}


gsl_vector *pget_min_distance(gsl_matrix *A, gsl_matrix *B, struct pind_func_dist *fparam)
{
    gsl_vector *acum = gsl_vector_alloc(A->size1);
    gsl_vector *euc = gsl_vector_alloc(B->size1);

    for (size_t s = 0; s < A->size1; s++) {

        gsl_matrix *Sc = repeatRowToSizeB(s,A,B);

        for (size_t s = 0; s < Sc->size1; s++) {
            gsl_vector_view cr = gsl_matrix_row(Sc,s);
            gsl_vector_view fcr = gsl_matrix_row(B,s);
            gsl_vector_set(euc,s,fparam->f(&cr.vector, &fcr.vector, fparam->pval));
        }

        gsl_vector_set(acum,s,gsl_vector_min(euc));
        gsl_matrix_free(Sc);
    }
    gsl_vector_free(euc);

    return acum;
}


double euclidian(gsl_vector *A, gsl_vector *B, double p)
{
    double res;
    gsl_vector *C = gsl_vector_alloc(A->size);
    gsl_vector_memcpy(C, A);
    gsl_vector_sub(C, B);

//    res = gsl_blas_dnrm2(C);
    gsl_vector_map(C,pin_pow,&p);
    res = pow(gsl_vector_sum(C),1.0/p);
    
    gsl_vector_free(C);
    return res;
}

double pin_dist_p(gsl_vector *A, gsl_vector *B, double p)
{
    double res;
    gsl_vector *C = gsl_vector_alloc(A->size);
    gsl_vector_memcpy(C, A);
    gsl_vector_sub(C, B);

    gsl_vector_map(C,pin_pow,&p);
    res = gsl_vector_sum(C);

    gsl_vector_free(C);
    return res;
}

double euclidianMax(gsl_vector *Z, gsl_vector *A, double inv)
{
    // m_p = p;
    // A - Z  IGD+
    // Z - A GD+
    double res;
    gsl_vector *C = gsl_vector_alloc(Z->size);
    gsl_vector_memcpy(C, (inv != 0)? A: Z);
    gsl_vector_sub(C, (inv != 0)? Z: A);

    for (size_t i = 0; i < C->size; i++) {
        double cval = gsl_vector_get(C,i);
        gsl_vector_set(C,i,(cval > 0)? cval : 0);
    }
    
    res = gsl_blas_dnrm2(C);
    gsl_vector_free(C);
    return res;
}

double GD(gsl_matrix *FP, gsl_matrix *S, double p)
{
    m_p = p;
    double indval = -1;

    struct pind_func_dist func = {euclidian, p};
    gsl_vector *acum = pget_min_distance(S,FP, &func);

    gsl_vector_map(acum, pin_pow, &p);
    indval = pow(gsl_vector_sum(acum), 1.0/p) / acum->size;

//    indval = gsl_vector_sum(acum) / acum->size;

    gsl_vector_free(acum);
    return indval;
}

double GDp(gsl_matrix *FP, gsl_matrix *S, double p)
{
    m_p = p;
    // foreach row
    double gdval = -1;

    struct pind_func_dist func = {euclidian, p};
    gsl_vector *acum = pget_min_distance(S,FP, &func);
    gdval = gsl_vector_sum(acum) / acum->size;

//    gsl_vector_map(acum,pin_pow,&p);

    gsl_vector_free(acum);
    return pow(gdval, 1.0/p);
}

double GDplus(gsl_matrix *FP, gsl_matrix *S, double p)
{
    m_p = p;
    double gdval = -1;

    struct pind_func_dist func = {euclidianMax, 0};
    gsl_vector *acum = pget_min_distance(S,FP, &func);

    gdval = gsl_vector_sum(acum) / acum->size;
    gsl_vector_free(acum);

    return gdval;
}

double IGD(gsl_matrix *FP, gsl_matrix *S, double p)
{
    double indval = -1;

    struct pind_func_dist func = {euclidian, p};
    gsl_vector *acum = pget_min_distance(FP,S, &func);

    gsl_vector_map(acum, pin_pow, &p);
    indval = pow(gsl_vector_sum(acum), 1.0/p) / acum->size;

    gsl_vector_free(acum);
    return indval;
}

double IGDp(gsl_matrix *FP, gsl_matrix *S, double p)
{
    m_p = p;
    double gdval = -1;

    struct pind_func_dist func = {euclidian, p};
    gsl_vector *acum = pget_min_distance(FP,S, &func);
    gdval = gsl_vector_sum(acum) / acum->size;

    gsl_vector_free(acum);
    return pow(gdval, 1.0/p);
}

// igd and igdplus are same
// make dummy function with funtion pointer
// to change operator
double IGDplus(gsl_matrix *FP, gsl_matrix *S, double p)
{
    m_p = p;
    double gdval = -1;

    struct pind_func_dist func = {euclidianMax, 1}; //pval controls inverse stat
    gsl_vector *acum = pget_min_distance(FP,S, &func);

    gdval = gsl_vector_sum(acum) / acum->size;
    gsl_vector_free(acum);

    return gdval;
}

double deltap(gsl_matrix *FP, gsl_matrix *S, double p)
{
    double gd = GDp(FP,S,p);
    double igd = IGDp(FP,S,p);
    return (gd > igd)? gd : igd;
}

int vector_dominate(gsl_vector *u, gsl_vector *v)
{
    // default u not dom v
    int all = 0;
    int any = 0;

    size_t size = u->size;
    double *ud = u->data;
    double *vd = v->data;
    //  return all(u .<= v,dims=2) .& any(u .< v, dims=2)

    for (size_t i = 0; i < size; ++i) {
        all += ud[i] > vd[i];
    }

    if (all == 0) {
        for (size_t i = 0; i < size; ++i) {
            any = (ud[i] < vd[i]) ? 1 : 0;
            if (any) break;
        }
    }

    return (any && all == 0)? -1 : 1;
}

bool uDomv(gsl_vector *u, gsl_vector *v)
{
    gsl_vector *w = gsl_vector_alloc(u->size);
    gsl_vector_memcpy(w, u);
    gsl_vector_sub(w, v);

    // all u <= v
    int all = 0;
    int any = 0;
    for (size_t i = 0; i < w->size; i++) {
        double cv = gsl_vector_get(w,i);
        all += !(cv <= 0)? 1: 0;
        any += (cv < 0)? 1: 0;
    }
    gsl_vector_free(w);
    return (all == 0 & any > 0);
}

double inSetTwoCover(gsl_matrix *B, gsl_matrix *A)
{
    double gdval = -1;
    double acum = 0;
    for (size_t b = 0; b < B->size1; b++) {

        gsl_matrix *Bc = repeatRowToSizeB(b,B,A);
        gsl_vector *euc = gsl_vector_alloc(Bc->size1);

        for (size_t b = 0; b < Bc->size1; b++) {
            gsl_vector_view br = gsl_matrix_row(Bc,b);
            gsl_vector_view ar = gsl_matrix_row(A,b);
            gsl_vector_set(euc,b,
                (vector_dominate(&ar.vector,&br.vector) == -1)? 0:1
                );
        }
        // printf("udomv %lf \n", gsl_vector_sum(euc));
        acum += (gsl_vector_sum(euc) > 0)? 1 : 0;
        // printf("acum %lf--\n", acum);

        gsl_vector_free(euc);
        gsl_matrix_free(Bc);
    }
    gdval = acum / B->size1;

    return gdval;
}

gsl_vector* inParetoRank(gsl_matrix *A, gsl_matrix *B)
{
    gsl_vector *rank = gsl_vector_alloc(A->size1);
    for (size_t a = 0; a < A->size1; a++) {

        gsl_matrix *Ac = repeatRowToSizeB(a,A,B);
        gsl_vector *euc = gsl_vector_alloc(Ac->size1);

        for (size_t ac = 0; ac < Ac->size1; ac++) {
            gsl_vector_view ar = gsl_matrix_row(Ac,ac);
            gsl_vector_view br = gsl_matrix_row(B,ac);
            gsl_vector_set(euc,ac,(vector_dominate(&ar.vector,&br.vector) == -1)? 0:1);
        }
        // it may be better to use a simple vector
        gsl_vector_add(rank,euc);
        gsl_vector_free(euc);
        gsl_matrix_free(Ac);
    }

    return rank;
}
