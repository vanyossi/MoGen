//
// Created by Iván Yossi on 30/04/22.
//

#include "mgn_initializer.h"

#include <stdbool.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>

#include "mgn_types.h"
#include "mgn_random.h"
#include "individual.h"

void mgn_init_rand(void* in, void* params)
{
    struct mgnp_init_params *pr = (struct mgnp_init_params*) params;
    mgnLimit* lim = pr->limit;

    size_t sizex = pr->ops->get_iparams(in).x->size;
    for (size_t i = 0; i < sizex; i++) {
        gsl_vector *x = pr->ops->get_iparams(in).x;
        gsl_vector_set(x, i, rnd_getUniform_limit(lim->min[i], lim->max[i]));
    }
}



// =========== Latin HyperCube ======

//static gsl_matrix *m_LHCm;
//static bool m_LHCinit = false;
//static size_t m_LHCcur;

struct pmgn_lhc_data_priv {
    size_t p_cur;
};

mgn_lhci*
mgn_init_new_lhci(size_t psize, size_t dim, mgnLimit *lim)
{
    mgn_lhci *lhci = calloc(1, sizeof(*lhci));
    lhci->psize = psize;
    lhci->dim = dim;
    lhci->limits = lim;
    lhci->m_points = gsl_matrix_alloc(psize,dim);

    lhci->_p = malloc(sizeof(lhci->_p));
    lhci->_p->p_cur = 0;

//    gsl_permutation *perm_row = gsl_permutation_calloc(psize);
//    gsl_ran_shuffle (rnd_get_generator(), perm_row->data, psize, sizeof(size_t));
//
//    for (size_t i = 0; i < lhci->m_points->size2; ++i) {
//        // each cuadrant
//        double delta = lim->max[i] / psize;
//        double spill = delta / 5;
//        double kmin = 0;
//
//        for (size_t j = 0; j < lhci->m_points->size1; ++j) {
//            double center = kmin + (lim->min[i] + delta) / 2;
//            gsl_matrix_set(lhci->m_points,j,i,
//                rnd_getUniform_limit(
//                    center - spill
//                    ,center + spill)
//            );
//            kmin += delta;
//        }
//        gsl_vector_view cvv = gsl_matrix_column(lhci->m_points,i);
//        gsl_permute_vector(perm_row, &cvv.vector);
//        gsl_ran_shuffle (rnd_get_generator(), perm_row->data, perm_row->size, sizeof(size_t));
//    }
//
//    gsl_permutation_free(perm_row);
    mgn_init_lhc_to_matrix(lhci->m_points, lim);

    return lhci;
}

void
mgn_init_lhc_to_matrix(gsl_matrix *m_a, mgnLimit *lim)
{
    size_t psize = m_a->size1;
    gsl_permutation *perm_row = gsl_permutation_calloc(psize);
    gsl_ran_shuffle (rnd_get_generator(), perm_row->data, psize, sizeof(size_t));

    for (size_t i = 0; i < m_a->size2; ++i) {
        // each cuadrant
        double delta = lim->max[i] / (double)psize;
        double spill = delta / 5;
        double kmin = 0;

        for (size_t j = 0; j < m_a->size1; ++j) {
            double center = kmin + (lim->min[i] + delta) / 2;
            gsl_matrix_set(m_a,j,i,
                rnd_getUniform_limit(
                    center - spill
                    ,center + spill)
            );
            kmin += delta;
        }
        gsl_vector_view cvv = gsl_matrix_column(m_a,i);
        gsl_permute_vector(perm_row, &cvv.vector);
        gsl_ran_shuffle (rnd_get_generator(), perm_row->data, perm_row->size, sizeof(size_t));
    }

    gsl_permutation_free(perm_row);
}

void mgn_init_lhc(void *in, void *lhcip)
{
    mgn_indv *indv = (mgn_indv*) in;
    mgn_indv_ops* ops = indv->ops;
    mgn_lhci *lhci = (mgn_lhci*) lhcip;

    gsl_vector *x = ops->get_iparams(in).x;

    gsl_matrix_get_row(x,lhci->m_points,lhci->_p->p_cur);
    lhci->_p->p_cur++;
}

void
mgn_lhci_reset(mgn_lhci *lhc)
{
    lhc->_p->p_cur = 0;
}

void
mgn_lhci_free(mgn_lhci *lhci)
{
    gsl_matrix_free(lhci->m_points);
    free(lhci->_p);
    free(lhci);
}


//void mgn_init_LHC_init(size_t psize, size_t dim, mgnLimit *lim)
//{
//
////    1st lower,upper/psize
////    2nd psize/psize, 2lower/psize
//    gsl_permutation *perm_row = gsl_permutation_calloc(psize);
//    gsl_ran_shuffle (rnd_get_generator(), perm_row->data, psize, sizeof(size_t));
//    m_LHCm = gsl_matrix_alloc(dim,psize);
//
//    for (size_t i = 0; i < m_LHCm->size1; ++i) {
//        // each cuadrant
//        double delta = lim->max[i] / psize;
//        double spill = delta / 5;
//        double kmin = 0;
//
//        for (size_t j = 0; j < m_LHCm->size2; ++j) {
//            double center = kmin + (lim->min[i] + delta) / 2;
//            gsl_matrix_set(m_LHCm,i,j,
//            rnd_getUniform_limit(
//                center - spill
//               ,center + spill)
//            );
////            printf("limit: %d, %d :: %.4f, %.4f %.6f\n",i, j, kmin + lim->min[i], kmin + delta, center - delta/4);
//            kmin += delta;
//        }
//        gsl_vector_view cvv = gsl_matrix_row(m_LHCm,i);
////        gsl_vector_fprintf(stdout, &cvv.vector, "%g");
//        gsl_permute_vector(perm_row, &cvv.vector);
////        gsl_ran_shuffle(rnd_get_generator(),cvv.vector.data,cvv.vector.size, sizeof(*cvv.vector.data));
//        gsl_ran_shuffle (rnd_get_generator(), perm_row->data, perm_row->size, sizeof(size_t));
//    }
//    m_LHCcur = 0;
//    m_LHCinit = true;
//}
//
//gsl_matrix* mgn_LHC_get()
//{
//    return m_LHCm;
//}

//void mgn_init_lhc(void *in, void* pop_iops)
//{
//    struct _mgn_i_ops *ops = (struct _mgn_i_ops*)pop_iops;
//    gsl_vector *x = ops->get_iparams(in).x;
//
//    gsl_matrix_get_col(x,m_LHCm,m_LHCcur);
//    m_LHCcur++;
//}

