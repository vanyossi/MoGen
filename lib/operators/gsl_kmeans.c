/*
         kmeans - K-Means Clustering Library for the GNU PSPP (People Should
Prefer PSPP) project.
         Copyright (C) 2011  Dr.Mehmet Hakan Satman <address@hidden>

         This program is free software: you can redistribute it and/or modify
         it under the terms of the GNU General Public License as published by
         the Free Software Foundation, either version 3 of the License, or
         (at your option) any later version.

         This program is distributed in the hope that it will be useful,
         but WITHOUT ANY WARRANTY; without even the implied warranty of
         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
         GNU General Public License for more details.

         You should have received a copy of the GNU General Public License
         along with this program.  If not, see <http://www.gnu.org/licenses/>.
     */


#include "gsl_kmeans.h"

struct Kmeans*
kmeans_create(double* data,
    int n,
    int m,
    int ngroups,
    int maxiter){
    int i,j;
    struct Kmeans *k = (struct Kmeans*) malloc(sizeof(struct Kmeans));
    k->data=gsl_matrix_alloc(n,m);
    k->centers=gsl_matrix_alloc(ngroups, m);
    k->ngroups=ngroups;
    k->index=gsl_vector_int_alloc(n);
    k->n=n;
    k->m=m;
    k->maxiter=maxiter;
    k->lastiter=0;
    for (i=0;i<n;i++){
        for(j=0;j<m;j++){
            gsl_matrix_set(k->data, i, j, data[i * m +j]);
        }
    }
    k->weights = (double*)malloc(sizeof(double) * k->index->size);
    k->v1 = gsl_vector_alloc(k->centers->size2);
    k->v2 = gsl_vector_alloc(k->centers->size2);
    k->v3 = gsl_vector_alloc(k->n);
    return(k);
}

void
kmeans_randomize_centers(struct Kmeans *kmeans){
    int i,j;
    double min=0,max=0;
    min=gsl_matrix_min(kmeans->data);
    max=gsl_matrix_max(kmeans->data);
    gsl_matrix_minmax(kmeans->data, &min, &max);
    for (i=0;i<kmeans->centers->size1;i++){
        for (j=0;j<kmeans->centers->size2; j++){
            gsl_matrix_set(kmeans->centers, i, j, min +
                                                  (((double)rand())/RAND_MAX)*(max-min));
        }
    }
}



void
kmeans_print(struct Kmeans* kmeans){
    int i,j;
    printf("Number of groups: %d\n",kmeans->ngroups);
    printf("Centers:\n");
    for (i=0;i<kmeans->centers->size1;i++) {
        for (j=0;j<kmeans->centers->size2;j++){
            printf("%f ",gsl_matrix_get(kmeans->centers, i,j));
        }
        printf("\n");
    }

    printf("Index:\n");
    for (i=0;i<kmeans->n;i++){
        printf("%d ",gsl_vector_int_get(kmeans->index, i));
    }
    printf("\nLast iter: %d\n",kmeans->lastiter);
}


void print_matrix(gsl_matrix *m){
    int i,j;
    for (i=0;i<m->size1;i++){
        for (j=0;j<m->size2;j++){
            printf("%f ",m->data[i * m->size2 + j]);
        }
        printf("\n");
    }
}


double
kmeans_euclidean_distance(gsl_vector *v1,
    gsl_vector *v2){
    double result=0.0;
    double val;
    int i;
    for (i=0;i<v1->size;i++){
        val=v1->data[i] - v2->data[i];
        result+=val*val;
    }
    return(result);
}



int
kmeans_num_elements_group(struct Kmeans *kmeans, int group){
    int total=0;
    int i;
    for (i=0;i<kmeans->n;i++){
        if(gsl_vector_int_get(kmeans->index,i)==group) total++;
    }
    return(total);
}


void
kmeans_recalculate_centers(struct Kmeans *kmeans){
    int i,j,h;
    int elm;
    double mean;
    gsl_vector *v1=kmeans->v3;

    for (i=0;i<kmeans->ngroups;i++){
        elm=kmeans_num_elements_group(kmeans,i);
        for (j=0;j<kmeans->index->size;j++){
            if(gsl_vector_int_get(kmeans->index,j)==i){
                kmeans->weights[j]=1.0;
            }else{
                kmeans->weights[j]=0.0;
            }
        }

        for (h=0;h<kmeans->m;h++){
            gsl_matrix_get_col(v1,kmeans->data, h);
            mean=gsl_stats_wmean(kmeans->weights, 1, v1->data ,1,
                v1->size);
            gsl_matrix_set(kmeans->centers, i,h, mean);
        }
    }
}


void
kmeans_calculate_indexes(struct Kmeans *kmeans){
    double dist;
    double mindist;
    int bestindex=0;
    int i,j;
    gsl_vector *v1 = kmeans->v1;
    gsl_vector *v2 = kmeans->v2;

    for (i=0;i<kmeans->n; i++){
        mindist=INFINITY;
        gsl_matrix_get_row(v1, kmeans->data, i);
        for (j=0;j<kmeans->ngroups; j++){
            gsl_matrix_get_row(v2, kmeans->centers,j);
            dist=kmeans_euclidean_distance(v1,v2);
            if(dist<mindist){
                mindist=dist;
                bestindex=j;
            }
        }
        gsl_vector_int_set(kmeans->index,i,bestindex);
    }
}



int
kmeans_check_converge(gsl_vector_int *current,
    gsl_vector_int *old){
    int i;
    int total=0;
    for (i=0;i<current->size;i++) {
        if(current->data[i] == old->data[i]) total++;
        old->data[i]=current->data[i];
    }
    return(current->size-total);
}



gsl_matrix*
kmeans_getGroup(struct Kmeans *kmeans, int index){
    int i;
    int j=0;
    int elm=kmeans_num_elements_group(kmeans,index);
    gsl_matrix *agroup=gsl_matrix_alloc(elm, kmeans->data->size2);
    gsl_vector *v1=gsl_vector_alloc(kmeans->data->size2);
    for(i=0;i<kmeans->data->size1;i++){
        if(kmeans->index->data[i]==index){
            gsl_matrix_get_row(v1, kmeans->data, i);
            gsl_matrix_set_row(agroup, j, v1);
            j++;
        }
    }
    return(agroup);
}



void
kmeans_cluster(struct Kmeans *kmeans){
    int i;
    gsl_vector_int *oldindex = gsl_vector_int_alloc(kmeans->index->size);
    kmeans_randomize_centers(kmeans);

    for (i=0;i<kmeans->maxiter;i++){
        kmeans->lastiter=i;
        kmeans_calculate_indexes(kmeans);
        kmeans_recalculate_centers(kmeans);
        if(kmeans_check_converge(kmeans->index, oldindex)==0) break;
    }

}
