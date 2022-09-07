#ifndef _MGN_RANDOM_
#define _MGN_RANDOM_

#include <gsl/gsl_rng.h>


void rnd_initialize();

void rnd_gen_free();

void rnd_set_seed(unsigned long seed);

unsigned long rnd_getUniform_int(int max);

double rnd_getUniform();

double rnd_getUniform_limit(double lower, double upper);

gsl_rng* rnd_get_generator();

#endif // _MGN_RANDOM_
