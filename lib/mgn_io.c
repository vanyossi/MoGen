/*
 *
 *  SPDX-FileCopyrightText: 2022$ Iván Yossi <ghevan@gmail.com>
 *
 *  SPDX-License-Identifier: GPL-2.0-or-later
 */

#include "mgn_io.h"


#include <stdio.h>
#include <stdlib.h>
//#include <string.h>

#define MAX_LINE_LENGTH 60

bool it_read_data(const char *filename, struct _inGroup_list* parsedData)
{
    FILE *istream;
    bool input_empty = false;
    char *strLine = (char*)calloc(sizeof(char), MAX_LINE_LENGTH);
    int linePointer;

    int res;
    double ndata;
    // printf("fname %s\n", filename);

    if (filename == NULL) {
        istream = stdin;
        // This wont work for piped streams
//        if (fseek(stdin, 0, SEEK_END), ftell(stdin) > 0)
//        {
//            rewind(stdin);
//        } else {
//            input_empty = true;
//        }

    } else {
        istream = fopen(filename, "r");
        input_empty = (istream)? false : true;
    }

    if(!input_empty) {

        // while(fgets(strLine, MAX_LINE_LENGTH, istream)) {
        // skip any comment starting line;
        char c;
        double num;
        int totalstep = 0;
        int step = 0;
        int stepi = 0;

        bool startCount = false;
        bool firstRun = true;
        bool readData = false;

        int col = 0;
        int row = 0;
        struct _inData_list* valueList = calloc(sizeof(struct _inData_list), 1);
        // if (scanf(strLine, "%[#/]", &c) == 1) {
        //     continue;
        // }

        // This does not work for stdin from pipes
        // better to read a buffer until newline character

//        while (!input_empty && fscanf(istream, "%[^\n] ", strLine) > EOF) {
        size_t bsize = MAX_LINE_LENGTH;
        ssize_t readsize = 0;
        while (readsize = getline(&strLine, &bsize,istream), readsize > 1) {
//            printf("size %zd", readsize);
            if (strLine[0] == '#') {
                if (readData) {
                    col = 0;
                    row = 0;
                    // save previous data and make a new one
                    inGroupList_insert(parsedData, valueList);
                    valueList = calloc(sizeof(struct _inData_list), 1);
                    readData = false; // if more than one # line
                }
                continue;
            }
//             printf("word %s", strLine);
            // startCount = (true & firstRun);
            row++;
            totalstep = 0;
            step = 0;
            while(sscanf(strLine + totalstep, "%lf%n", &num, &step) > 0) {
                col++;
                totalstep += step;
                // printf("num %lf::%d %d\n", num, step, totalstep);
                inData_insert_value(valueList, num);
//                printf("value %.4f, ", num);
            }
//            printf("\n");
            // printf("num %lf::%ld\n", strtof(strLine, NULL), ftell(istream));
            // totalstep += step;
            // printf("-- %d %d\n", totalstep, step);
            // stepi++;
            // if (stepi >= 100) {
            //     break;
            // }
            valueList->col = col / row;
            valueList->row = row;
            readData = true;
        }
        // insert last datagroup
        inGroupList_insert(parsedData, valueList);
        // printf("col, row (%d,%d)\n", col / row, row);
        // }
        // res = fscanf(istream, "%lf", &ndata);

        fclose(istream);
    }
    free(strLine);
    return !input_empty;
}
