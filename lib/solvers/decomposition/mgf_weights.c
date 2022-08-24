/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

// TODO 1: migrate to mogen
// TODO 2: design a better error handling strat
//      maybe global with static struct for report (?)
#include "mgf_weights.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mgf_global_types.h"

/**
 * Set the parameters for the Simplex-Lattice Design in a layer configuration approach
 * @param wp       The structure with the layers configuration
 * @param layers_size   The number of layer
 * @param ...      The value for each layer
 * @return         The total number of layers
 */
int WS_SLD_config_instance(WeightParams *wp, int layers_size, int layers[layers_size])
{
    wp->layers = layers_size;
    wp->H_layer = (int*) malloc(sizeof(int) * (unsigned int)layers_size);

    for (int i = 0; i < layers_size; ++i)
    {
        wp->H_layer[i] = layers[i];
        //printf("H%d: %d\n", i + 1, wp.H_layer[i]);
    }
    return wp->layers;
}

/**
 * Simplex-Lattice Design layers configuration (as NSGA3)
 * @param ws    Weight set parameters
 * @return      The size of the population according to the configuration
 */
unsigned int WS_SLD_config(WeightParams *wp, int nobj)
{
    unsigned int psize = 0;
    int L;

    L = (nobj == 2) ? WS_SLD_config_instance(wp, 1, (int[]){100}) : /* 2objs (1 layer): 99=>N=100 */
        (nobj == 3) ? WS_SLD_config_instance(wp, 1, (int[]){19}) : /* 3objs (1 layer): 19=>N=210 */
        (nobj == 4) ? WS_SLD_config_instance(wp, 1, (int[]){9}) : /* 4objs (1 layer): 9=>N=220 */
        (nobj == 5) ? WS_SLD_config_instance(wp, 1, (int[]){6}) : /* 5objs (1 layer): 6=>N=210 */
        (nobj == 6)? WS_SLD_config_instance(wp, 2, (int[]){4, 3}) : /* 6objs (2 layer): 4=>126 + 3=>56 => N=182 */
        //(nobj == 6) ? EMOA_WS_SLD_cfg_instance(ws, 1, 7) : /* 6objs (2 layer): 4=>126 + 3=>56 => N=182 */
        (nobj == 7) ? WS_SLD_config_instance(wp, 2, (int[]){4, 2}) : /* 7objs (2 layer): 4=>210 + 2=>28 => N=238 */
        (nobj == 8) ? WS_SLD_config_instance(wp, 2, (int[]){3, 2}) : /* 8objs (2 layer): 3=>120 + 2=>36 => N=156 */
        (nobj == 9) ? WS_SLD_config_instance(wp, 2, (int[]){3, 2}) : /* 9objs (2 layer): 3=>165 + 2=>45 => N=210 */
        (nobj == 10) ? WS_SLD_config_instance(wp, 2, (int[]){3, 2}) : /* 10objs (2 layer): 3=>220 + 2=>55 => N=275 */
        (nobj == 11) ? WS_SLD_config_instance(wp, 2, (int[]){2, 2}) : /* 11objs (2 layer): 2=>66 + 2=>66 => N=132 */
        (nobj == 12) ? WS_SLD_config_instance(wp, 2, (int[]){2, 2}) : /* 12objs (2 layer): 2=>78 + 2=>78 => N=156 */
        (nobj == 13) ? WS_SLD_config_instance(wp, 2, (int[]){2, 2}) : /* 13objs (2 layer): 2=>91 + 2=>91 => N=182 */
        (nobj == 14) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 14objs (2 layer): 2=>105 + 1=>14 => N=119 */
        (nobj == 15) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 15objs (2 layer): 2=>120 + 1=>15 => N=135 */
        (nobj == 16) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 16objs (2 layer): 2=>136 + 1=>16 => N=152 */
        (nobj == 17) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 17objs (2 layer): 2=>153 + 1=>17 => N=170 */
        (nobj == 18) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 18objs (2 layer): 2=>171 + 1=>18 => N=189 */
        (nobj == 19) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 19objs (2 layer): 2=>190 + 1=>19 => N=209 */
        (nobj == 20) ? WS_SLD_config_instance(wp, 2, (int[]){2, 1}) : /* 20objs (2 layer): 2=>210 + 1=>20 => N=230 */
        -1;
    assert(L > 0);
    UNUSED(L);

    //printf("flag: %d\n", layers);
    //printf("ws.layers: %d\n", ws->layers);
    for (int i = 0; i < wp->layers; ++i)
    {
        psize += EMOA_combination(nobj, wp->H_layer[i]);
        //printf("H%d: %d\n", i + 1, ws->H_layer[i]);
    }
//    printf("N=%d\n", psize);
    assert(psize > 0);
    return psize;
}


/**
 * File configuration
 * @param wp    Weight set parameters
 * @param option The possible option to set:
 *          1 According to the simplex-lattice layer configuration
 *          2 Resolution 100*k (k: number of objective functions)
 *          3 Constant population = 210
 * @return The size of the population according to the configuration
 */
unsigned int WS_FILE_settings(WeightParams *wp, int nobj, WeightFileOption option)
{
    unsigned int psize = 0;

    wp->layers = 0;
    switch (option)
    {
        case W_SIMPLE_LATICE: // Population size according to the layers in simplex-lattice
            psize = (nobj == 2) ? 100 :
                (nobj == 3) ? 210 :
                (nobj == 4) ? 220 :
                (nobj == 5) ? 210 :
                (nobj == 6) ? 182 :
                (nobj == 7) ? 238 :
                (nobj == 8) ? 156 :
                (nobj == 9) ? 210 :
                (nobj == 10) ? 275 :
                (nobj == 11) ? 132 :
                (nobj == 12) ? 156 :
                (nobj == 13) ? 182 :
                (nobj == 14) ? 119 :
                (nobj == 15) ? 135 :
                (nobj == 16) ? 152 :
                (nobj == 17) ? 170 :
                (nobj == 18) ? 189 :
                (nobj == 19) ? 209 :
                (nobj == 20) ? 230 : 0;
            break; //
        case W_RES_100K: // Population size with resolution 100*k (k: number of objectives)
            psize = 100U * nobj;
            break;
        case W_RES_210K: // Population size with resolution 210 (k: number of objectives)
            psize = 210U * nobj;
            break;
        case W_DEFAULT:
            printf("ERROR: weight File option %d Option no supported\n", W_DEFAULT);
            exit(0);
    }
    assert(psize > 0);
    return psize;
}


/**
 * Set the weights
 * @param ws        Weight set parameters
 * @param foption   Weight generator for file option
 *                      SIMPLE_LATICE
 *                      RES_100K
 *                      RES_210K
 * @return          size of popuation
 */
unsigned int weights_settings(int nobj, WeightParams *wp, WeightFileOption foption)
{
    unsigned int psize = 0;

    psize = (wp->method == WM_WEIGHT_SLD) ?
            WS_SLD_config(wp, nobj) : WS_FILE_settings(wp, nobj, foption);

    sprintf(wp->file, "WEIGHTS/SOBOL-%dobjs-%dwei.ws", nobj, psize);

    assert(psize > 0);
    return psize;
}

/**
 * Recursive generation only one layer
 * @param W     The weight set
 * @param w     The weight vector
 * @param p_id  The id of the weight vector
 * @param h_top The top of parameter H
 * @param H     The parameter H for the simplex-lattice design
 * @param wj    The axis in the recursion
 */
static void WS_recursive_one_layer(struct mgf_weight_set *ws, double *w, unsigned int *p_id, int h_top, int H, int wj, int nobj)
{
    int i;
    if (wj == (nobj - 1)) {
        w[wj] = (1.0 * h_top) / (double) H;
        memcpy(ws->weight[*p_id].w, w, sizeof(double) * nobj);
        *p_id = *p_id + 1;

    } else {
        for (i = 0; i <= h_top; i++) {
            w[wj] = (1.0 * i) / (double) H;
            WS_recursive_one_layer(ws, w, p_id, (h_top - i), H, (wj + 1), nobj);
        }
    }
}

/**
 * Initialize Weight set layers
 * @param W The weight_set
 * @param ws Weight set parameters
 */
static void initialize_WS_SLD_Layers(WeightParams *wp, int nobj)
{
    int i;
    unsigned int j, l, l1;
    unsigned int p_id = 0;
    int H;

    struct mgf_weight_set *W = wp->wset;
    double *w = calloc(nobj, sizeof(double));

    double centroid = 1.0 / (double) nobj;

    for (i = 0; i < wp->layers; ++i)
    {
        H = wp->H_layer[i];
        if (i == 0)
        {
            assert(p_id < W->size);
            //set_double_array(w, 0.0, mop.nobj);
            for (l = 0; l < nobj; ++l)
            {
                w[l] = 0.0;
            }
            WS_recursive_one_layer(wp->wset, w, &p_id, H, H, 0, nobj);
        }
        else
        {
            l1 = p_id;
            assert(p_id < W->size);
            //set_double_array(w, 0.0, mop.nobj);
            for (l = 0; l < nobj; ++l)
            {
                w[l] = 0.0;
            }
            WS_recursive_one_layer(wp->wset, w, &p_id, H, H, 0, nobj);
            for (l = l1; l < p_id; ++l)
            {
                for (j = 0; j < nobj; ++j)
                {
                    wp->wset->weight[l].w[j] = (centroid + wp->wset->weight[l].w[j]) / 2.0;
                }
            }
        }
    }
    assert(W->size == p_id);
    UNUSED(W);

    /* Labeling layers */
    for (i = 0; i < (int)wp->wset->size; i++)
    {
        wp->wset->weight[i].id = i;
    }
    free(w);
    return;
}

/***
 * Initialize the weight vectors according from the file (file_method)
 * @param W   The weight set
 * @param ws  Weight set parameters
 */
static void initialize_WS_FILE(WeightParams *wp, int nobj)
{
    unsigned int i, j;
    char file_name[1024];
    FILE *W_file = NULL;

    strcpy(file_name, wp->file);
    W_file = (FILE*) fopen(file_name, "r");
    if (W_file == NULL)
    {
        printf("ERROR: The %s file does not exist (initialize_WS_FILE).\n", file_name);
        exit(0);
    }
    /* vprint("Charging weight vectors from file: %s\n", file_name); */

    for (i = 0; i < wp->wset->size; ++i)
    {
        for (j = 0; j < nobj; ++j)
        {
            // TODO: use strtod instead (?)
            if (fscanf(W_file, "%lf", &wp->wset->weight[i].w[j]) == EOF)
            {
                printf("ERROR: Cannot read %s (initialize_weights_UD).\n", file_name);
                exit(0);
            }
        }
        // Labeling layers
        wp->wset->weight[i].id = (int)i;

    }
    fclose(W_file);

    return;
}

/**
 * Save the weight set in a file
 * @param W        The weight set
 * @param fname    The file name
 */
static void WS_save(struct mgf_weight_set *ws, char *fname, int nobj)
{
    unsigned int i, j;
    FILE *W_file;

    W_file = fopen(fname, "w+");
    for (i = 0; i < ws->size; ++i) {
        for (j = 0; j < nobj; ++j) {
            fprintf(W_file, "%lf ", ws->weight[i].w[j]);
        }
        fprintf(W_file, "\n");
    }

    fclose(W_file);
    return;
}

/**
 * Initialize the weight set
 * @param W         The weight set
 * @param ws        Weight parameter options
 */
void weights_initialize_set(WeightParams *wp, int nobj) {
    if(wp->method == WM_WEIGHT_SLD){
        initialize_WS_SLD_Layers(wp, nobj);

    } else if (wp->method == WM_WEIGHT_FILE){
        initialize_WS_FILE(wp, nobj);
    }
    WS_save(wp->wset, "current_WS.dat", nobj);
    return;
}
