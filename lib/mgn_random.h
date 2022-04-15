#ifndef _MGN_RANDOM_
#define _MGN_RANDOM_

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
    // gsl_rng_set(r, seed)
    return;
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



#endif // _MGN_RANDOM_