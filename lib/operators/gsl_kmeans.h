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

#ifndef MOGEN_GSL_KMEANS_H
#define MOGEN_GSL_KMEANS_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

/*
Struct KMeans:
Holds all of the information for the functions.
*/
struct Kmeans {
    gsl_matrix *data;           //User Data (Given by the user)
    gsl_matrix *centers;        //Centers for groups
    gsl_vector_int *index;      //Integer values from zero to ngroups.
                                //    Shows group of an observation.
    gsl_vector *v1,*v2,*v3;     //Temporary vector for program. Do not use.
    int ngroups;                //Number of group. (Given by the user)
    int n;                      //Number of observations. (Given by theuser)
    int m;                      //Number of observations. (Given by theuser)
    int maxiter;                //Maximum number of iterations (Given bythe user)
    int lastiter;               //Show at which iteration it found thesolution.
    double *weights;            //Double values for handling weights forprogram use.
};

/*
Creates a Kmeans structure for a given data 'data', number of observations
'n',
number of columns 'm', number of groups 'ngroups' (it is usually called
'k') and
number of maximum iterations 'maxiter'.
*/
struct Kmeans* kmeans_create(double* data, int n, int m, int ngroups, int
maxiter);


/*
Randomly chooses centers of 'ngroup' groups.
These centers are initial and will be changed iteratively.
*/
void kmeans_randomize_centers(struct Kmeans *kmeans);

/*
Prints the given Kmeans structure
*/
void kmeans_print(struct Kmeans* kmeans);


/*
Calculates the squared euclidean distance between vector v1 and v2
*/
double kmeans_euclidean_distance(gsl_vector *v1, gsl_vector *v2);


/*
Calculates and returns the number of elements contained in the
specific group.
*/
int kmeans_num_elements_group(struct Kmeans *kmeans, int group);


/*
Calculates group centers. Those centers are calculated using iteration
and they are usually different from the initial centers. Recalculation
process remains while the last two solutions are not equal.
*/
void kmeans_recalculate_centers(struct Kmeans *kmeans);


/*
Constructs the index variable. This variable shows the current
group of the each single observation.
*/
void kmeans_calculate_indexes(struct Kmeans *kmeans);


/*
Checks if the last two index variables are equal. If they are equal,
algorithm can not find a better classification anymore and stops.
*/
int kmeans_check_converge(gsl_vector_int *current, gsl_vector_int *old);


/*
This is the main method of the algorithm.
*/
void kmeans_cluster(struct Kmeans *kmeans);


#endif //MOGEN_GSL_KMEANS_H
