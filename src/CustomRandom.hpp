//Implements initialization of GSL Mersenne Twister

#ifndef CUSTOMRANDOM_HPP
#define CUSTOMRANDOM_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class CustomRandom
{
public:
    //Returns an instance of GSL's default random number generator (rng), which implements the Mersenne Twister (MT19937), seeded with seed
    static gsl_rng* init_rng(long unsigned int seed);
};

#endif