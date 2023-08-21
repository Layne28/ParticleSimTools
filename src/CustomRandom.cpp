#include "CustomRandom.hpp"

gsl_rng* CustomRandom::init_rng(long unsigned int seed)
{
    gsl_rng *rg;
 	const gsl_rng_type *T;
	T = gsl_rng_default;
	rg = gsl_rng_alloc(T);
	srand((unsigned) seed);
	gsl_rng_env_setup();
	gsl_rng_set(rg, seed);
    return rg;
}