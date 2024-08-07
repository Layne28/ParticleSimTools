//A LabBench consists of (1) a System, (2) a Solver, (3) an Observer, and (4) a ParamDict.
//The LabBench uses the ParamDict to read in and store parameters.
//It creates a System from the parameters, and uses a Solver to advance the System in time.
//The Observer writes data from the System to file over time.
//A LabBench optionally takes in a filename and a seed.
//If the filename is not empty, then the ParamDict reads in parameters from the corresponding file.
//If the seed is specified, then LabBench uses that number to seed the RNG. Otherwise, it uses a seed of 1.

#ifndef LABBENCH_HPP
#define LABBENCH_HPP

#include <string>

#include "System.hpp"
#include "Solver.hpp"
#include "Observer.hpp"
#include "ParamDict.hpp"
#include "boost/graph/adjacency_list.hpp"

class LabBench
{
public:
    ParamDict params;
    System sys;
    Observer obs;
    Solver solver;
    gsl_rng *rg;

    std::string experiment = "standard";
    int randomization_steps = 0;
    int equil_steps = 0;
    int production_steps = 0;
    int info_freq = 1;
    int do_energy_minimize = 1;

    /*** Methods ***/

    //constructor
    LabBench(ParamDict& theParams, gsl_rng*& theGen);

    //destructor
    ~LabBench();

    //run trajectory
    void run(int nstps=-1, std::string subdir="/", int net_freq=-1, int therm_freq=-1, int noisegen_freq=-1);

    //ffs functions
    auto run_ffs_stage1(int N0, double l0, double la, double lb);

    //experiments
    void do_experiment(std::string expt); //TODO: keep this Public, move rest to private
    void run_standard_experiment(); //equilibrate and then observe an unperturbed system
    void run_ffs_experiment();

    void run_swollen_experiment();
    void run_compression_experiment();
    void run_shear_experiment();
    void run_sine_perturb_experiment();
};

#endif