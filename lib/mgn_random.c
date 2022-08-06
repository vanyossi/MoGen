//
// Created by Iván Yossi on 20/04/22.
//

#include "mgn_random.h"

#include <math.h>
#include <gsl/gsl_rng.h>

static gsl_rng *randGenerator;

void rnd_initialize()
{
    if (randGenerator) { gsl_rng_free(randGenerator); }

    const gsl_rng_type * randType;
    gsl_rng_env_setup();
    randType = gsl_rng_default;
    randGenerator = gsl_rng_alloc(randType);

    return;
}

void rnd_set_seed(unsigned long seed)
{
    gsl_rng_set(randGenerator, seed);
}

unsigned long rnd_getUniform_int(int max)
{
    if (!randGenerator) { rnd_initialize(); }
    return gsl_rng_uniform_int(randGenerator, max);
}

double rnd_getUniform()
{
    if (!randGenerator) { rnd_initialize(); }
    return gsl_rng_uniform(randGenerator);
}

double rnd_getUniform_limit(double lower, double upper)
{
    if (!randGenerator) { rnd_initialize(); }

    double value;
    double diff = fabs(upper - lower);

    value = gsl_rng_uniform(randGenerator) * diff + lower;

    return value;
}

gsl_rng* rnd_get_generator()
{
    if (!randGenerator) { rnd_initialize(); }
    return randGenerator;
}
