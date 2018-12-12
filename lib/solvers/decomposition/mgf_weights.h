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


#ifndef MOGEN_MGF_WEIGHTS_H
#define MOGEN_MGF_WEIGHTS_H

enum weight_method_t {
  WM_WEIGHT_SLD,
  WM_WEIGHT_FILE
};

typedef enum weight_foption {
    W_SIMPLE_LATICE,
    W_RES_100K,
    W_RES_210K,
    W_DEFAULT,
} WeightFileOption;

struct mgf_weight {
    int id;       //!< Id
    double *w;    //!< Weight
};

struct mgf_weight_set {
    int size;       //!< Id
    struct mgf_weight *weight;    //!< Weight
};

typedef struct
{
    int transformation;     //!< Any transformation to the weight set
    char file[1024];        //!< File name of weights set
    int layers;             //!< Number of Layers for Simplex-Lattice design
    enum weight_method_t method;             //!< Method for the weights
    int *H_layer;           //!< Layer value for the Simplex-Lattice design
    struct mgf_weight_set *wset;
} WeightParams;

// neighborhood
struct mgf_neighborhood {
    unsigned int size;
    int *members;
};

struct mgf_neighborhood_set {
    unsigned int size;
    struct mgf_neighborhood *neighborhood;
};


unsigned int weights_settings(int nobj, WeightParams *wp, WeightFileOption foption);

void weights_initialize_set(WeightParams *wp, int nobj);

//void moead_neighborhood_set(int nobj, struct mgf_neighborhood_set *ns, struct mgf_weight_set *ws);


//moead.weights.transformation = 0;
//moead.weights.method = (lmop->nobj < 20) ? WEIGHTS_SLD : WEIGHTS_FILE;
//moead.pop_size = weights_settings(&moead.weights, RES_210K);


#endif //MOGEN_MGF_WEIGHTS_H
