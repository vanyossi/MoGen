//
// Created by Iván Yossi on 29/07/22.
//

#ifndef MOGEN_MGN_TEST_H
#define MOGEN_MGN_TEST_H

#include <stdbool.h>
#include <stdio.h>

static int tests_run = 0;
static int tests_pass = 0;

bool mgn_test(char* name, bool (*func)())
{
    bool pass;
    tests_run++;

    printf("%s \t", name);
    if(func()){
        tests_pass++;
        pass = true;
        printf("PASS!\n");
    } else {
        pass = false;
        printf("FAIL!\n");
    }
    return pass;
}

#endif //MOGEN_MGN_TEST_H
