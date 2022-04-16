#ifndef _MGN_MOEAD_SOLVER_
#define _MGN_MOEAD_SOLVER_

#include <stdlib.h>
#include <stdbool.h>

#include "mgn_moa.h"
#include "population.h"
#include "individual.h"

struct moead_features {
    double some;
};

static struct moead_features feat;

void moead_run(mgnMop *mop)
{
    /* general framework assuming Tcheby
     *  sea W un conjunto de N vectores de peso
     *  sea zmin el punto de referencia
     *
     *  moead minimiza en una sola corrida N problemas
     *
     *  se define un vecindario para W_i como un conjunto de
     *      de los vectore de peso más cercanos
     *
     *  Población las mejores soluciones hasta ahora.
     *  En cada gen tenemos
     *      |P| = N
     *      FV_i = Fitnes values of x^i
     *      z = el mejor valor encontrado para f_i
     *      PE = población externa para almacenar soluciones no dom
     *
     *
     *  --ALGO
     *  input(mop, stop, N subprob, W, T = num vectors en vecidario)
     *  1.1 EP = 0
     *  1.2 Calcula euclidean entre los vectores y encuentra los T mas
     *      cercanos.
     *  1.3 Genera la población inicial t obten FV
     *  1.4 init z
     *
     *  2 Actualiza
     *  2.1 Reproducción, select random j,k y reproducirlos -> y
     *  2.2 Mejora:, aplica heuristica de mejora para y -> y^i
     *  2.3 Actualiza z, por cada j en m, para y^i z_j = f(y^i)
     *  2.4 Actualiza el vecindario.
     *  2.5 Actualiza EP
     *      remover todos las solcuiones dom por f(y^i)
     *      alade f(y^i) si nadie le domina.
     */

    /*
     * calculamos la distancia de todos los vectores contra todos (matrix)
     * Para W_i obtenemos los T más cercanos.
     *
     */

}

bool moead_stop()
{
    return false;
}

mgnMoa* moead_init(){
    mgnMoa* moead = (mgnMoa*)calloc(1, sizeof(mgnMoa));
    moead->name = "MOEA/D";
    moead->run = moead_run;
    moead->stop = moead_stop;

    return moead;
}



#endif // _MGN_MOEAD_SOLVER_
