//ActiveSolver is a child of the Solver class.
//It creates an active noise Generator and uses it
//to add spatiotemporally correlated random "active"
//forces to the network in a Brownian dynamics scheme.

#ifndef ACTIVESOLVER_HPP
#define ACTIVESOLVER_HPP

#include "Solver.hpp"
#include <fftw3.h>
#include <complex>
#include <string>
#include <sstream>
#include <angen/Generator.hpp>

class ActiveSolver : Solver
{
public:

    /*** Variables ***/
    double va = 1.0;
    double force_thresh = 200.0;

    /*** Methods ***/

    //constructor (needs to take System, unlike Solver)
    ActiveSolver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~ActiveSolver();

    //Take a step forward in time
    void update(System &theSys);
    void update_adaptive(System &theSys, double deet, int level);

    //use fft to compute real-space active noise on grid
    std::vector<arma::vec> get_active_noise_forces(System &theSys, Generator &gen);
    std::vector<arma::vec> get_thermal_forces(System &theSys, double deet);

private:
    Generator* anGen; //active noise generator
};

#endif