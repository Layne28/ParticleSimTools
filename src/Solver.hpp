//A Solver advances a System in time according to some (for now, configuration-space) dynamics. 
//For now, it only supports overdamped Langevin dynamics solved via the Euler-Maruyama method.

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <fftw3.h>
#include <complex>
#include <string>
#include <sstream>
#include "Solver.hpp"
#include "CustomRandom.hpp"
#include "System.hpp"
#include <angen/Generator.hpp>

class Solver
{
public:

    /*** Variables ***/
    double dt;    //timestep
    double gamma; //friction
    int do_subtract_com = 0; //If =1, subtract center of mass force
                             //from random forces
    int do_adaptive_timestep = 0;

    double force_thresh = 200.0; //adjust dt if force is higher than this

    //Active noise parameters
    double va = 1.0;
    
    gsl_rng *rg;

    /*** Methods ***/

    //constructor
    Solver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~Solver();

    //take a step forward in time
    void update(System &theSys);
    void update_adaptive(System &theSys, double deet, int level);

    std::vector<arma::vec> get_thermal_forces(System &theSys, double deet);
    std::vector<arma::vec> get_aoup_forces(System &theSys);
    std::vector<arma::vec> get_subtracted_com_forces(System &theSys, std::vector<arma::vec> forces);
    void update_self_prop_vel(System &theSys, int index, double deet);

    //active noise methods
    std::vector<arma::vec> get_active_noise_forces(System &theSys, Generator &gen);

private:
    Generator* anGen; //active noise generator
};

#endif