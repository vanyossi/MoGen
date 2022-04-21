#ifndef _MGN_RANDOM_
#define _MGN_RANDOM_

#include <gsl/gsl_rng.h>


void rnd_initialize();

unsigned long rnd_getUniform_int(int max);

double rnd_getUniform();

double rnd_getUniform_limit(double lower, double upper);

gsl_rng* rnd_get_generator();

#endif // _MGN_RANDOM_
