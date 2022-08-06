/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#ifndef MOGEN_MGN_IO_DATALIST_H
#define MOGEN_MGN_IO_DATALIST_H

#include <gsl/gsl_matrix.h>

typedef struct _inGroup_list mgn_io_datagroup;
typedef struct _inData_list mgn_io_data;

//TODO make some of this private

struct _inData
{
    double value;
    struct _inData* next;
};

struct _inData_list
{
    int col;
    int row;
    struct _inData* next;
};

struct _inGroup
{
    struct _inData_list* data;
    struct _inGroup* next;
};

struct _inGroup_list
{
    int size;
    struct _inGroup* next;
};

void inData_insert_value(struct _inData_list* list, double value);
void inData_delete(struct _inData_list* list);
void inData_toArray(struct _inData_list *list, double **array, int length);

void inGroupList_insert(struct _inGroup_list* group, struct _inData_list* list);
void inGroupList_delete(struct _inGroup_list* group);

void inData_printValues(struct _inData_list *list);
void inGroup_printData(struct _inGroup_list *group);

struct _inData_list* inGroup_getListAt(struct _inGroup_list *group, int index);


// GSL data transform types
gsl_matrix* inData_toGSLMatrix(struct _inData_list* list);


#endif //MOGEN_MGN_IO_DATALIST_H
