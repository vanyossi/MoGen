/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_io_datalist.h"

#include <stdlib.h>
#include <stdio.h>

void inData_insert_value(struct _inData_list* list, double value)
{
    struct _inData *new = malloc(sizeof(struct _inData));
    new->value = value;
    new->next = 0;

    if (list->next != 0) {
        struct _inData *current = list->next;
        while (current->next != 0) {
            current = current->next;
        }
        current->next = new;
    } else {
        list->next = new;
    }
    return;
}

void inData_delete(struct _inData_list* list)
{
    struct _inData *prev = 0;
    struct _inData *current = list->next;
    while (current != 0) {
        prev = current;
        current = current->next;
        free(prev);
    }
    return;
}

void inData_toArray(struct _inData_list *list, double **array, int length)
{
    double *data = *array;
    struct _inData *current = list->next;

    int i = 0;
    while (current != 0 & i < length) {
        data[i++] = current->value;
        current = current->next;
    }

    return;
}

gsl_matrix* inData_toGSLMatrix(struct _inData_list *list)
{
    int row = list->row;
    int col = list->col;
    gsl_matrix *m = gsl_matrix_alloc(row, col);

    struct _inData *current = list->next;

    for (size_t i = 0; i < row; i++) {
        for (size_t k = 0; k < col; k++) {
            gsl_matrix_set(m,i,k,current->value);
            current = current->next;
        }
    }
    return m;
}

void inGroupList_insert(struct _inGroup_list* group, struct _inData_list* list)
{
    struct _inGroup *new = malloc(sizeof(struct _inGroup));
    new->data = list;
    new->next = 0;

    struct _inGroup *current = group->next;
    if (group->size != 0) {
        while (current->next != 0) {
            current = current->next;
        }
        current->next = new;
    } else {
        group->next = new;
    }
    group->size++;

    return;
}

void inGroupList_delete(struct _inGroup_list* group)
{
//    void* current = pop->first;
//    mgnt_pop_free(ind_free) = pop->ops->free;
//    if(current != 0){
//        struct _mgn_i_ops *ops = pop->ops->get_iops(current);
//        while (ops->next(current) != 0) {
//            void* for_free = current;
//            current = ops->next(current);
//            ind_free(for_free);
//            free(for_free);
//        }
//        ind_free(current);
//        free(current);
//    }
//    ind_free(pop->I);
//    free(pop->I);
//    free(pop);

    struct _inGroup *prev = 0;
    struct _inGroup *current = group->next;
    while (current != 0) {
        inData_delete(current->data);
        free(current->data);
        prev = current;
        current = current->next;

        free(prev);

//        while (current->next != 0) {
//            prev = current;
//            current = current->next;
//            inData_delete(prev->data);
//            free(prev);
//        }

    }
    return;
}



void inData_printValues(struct _inData_list *list)
{
    struct _inData *current = list->next;
    int col = list->col;
    int ccol = 0;
    while (current != 0) {
        printf("%lf, ", current->value);
        ccol++;

        if (col == ccol) {
            printf("\b\b \n");
            ccol= 0;
        }
        current = current->next;
    }
    printf("\b\b \n");
    return;
}

void inGroup_printData(struct _inGroup_list *group)
{
    struct _inGroup *current = group->next;
    while (current != 0) {
        inData_printValues(current->data);
        current = current->next;
    }
    return;
}

struct _inData_list* inGroup_getListAt(struct _inGroup_list *group, int index)
{
    struct _inData_list* res = 0;

    int pos = 0;
    if (index <= group->size) {
        struct _inGroup *current = group->next;
        while (pos < index) {
            current = current->next;
            pos++;
        }
        res = current->data;
    }
    return res;
}
