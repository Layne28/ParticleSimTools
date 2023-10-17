//A System consists of a set of Particles along with boundary conditions and control parameters (temperature, external field, etc.)

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include "Particle.hpp"
#include "ParamDict.hpp"
#include "Observer.hpp"
#include "CustomRandom.hpp"
#include "../../Common-Tools/neighborgrid.h"

class Observer;

class System
{
private:
    Observer *obs; //use to write output
    gsl_rng *rg;

public:
    //For now we assume NVT ensemble
    int N;      //No. of particles
    double rho; //density
    double a;   //lattice constant
    double phi = -1.0; //packing fraction
    double kT;  //temperature
    int dim;    //no. of spatial dimensions
    double dt;  //timestep

    int is_network = 0;     //Whether particles are connected to other particles by springs
    int can_bonds_break = 0; //Whether network connectivity can change
    int is_aoup = 0; //Whether to make particles AOUPs

    double rcut; //cutoff distance for pair potential
    double sigma;
    double epsilon;
    double l0; //spring rest length
    double K; //spring constant
    double aoup_D0; //aoup diffusion constant
    double aoup_tau; //aoup correlation time
    double drmax;
    std::string bonded_potential_type = "harmonic";
    std::string nonbonded_potential_type = "none";
    std::string particle_protocol;
    std::string spring_protocol;

    int time; //No. of timesteps taken (can be reset)

    std::vector<Particle> particles;     //all the particles in the system
    std::vector<std::vector<int>> image; //which periodic image each particle is in (starts at (0,0))
    std::vector<double> edges;           //box dimensions
    std::vector<int> unit_cells;         //no. of unit cells in each direction (for lattice)
    std::vector<int> is_periodic;        //periodic or not in each dimension

    int do_cell_list = 1;     //whether to use the cell list for computing forces
    int do_neighbor_grid = 0; //whether to use the neighbor grid for computing forces

    //Cell list data structures (assume 2d for now)
    //See A&T 5.3.2
    int ncell_x, ncell_y;
    double cellsize_x, cellsize_y;
    int *head;       //first particle in each cell
    int *list;       //list of particles in cell
    int *cellndx;    //which cell each particle belongs to
    int **cellneigh; //each cell keeps list of its neighboring cells


    NeighborGrid<Particle, 2> *grid;

    /*** Methods ***/

    //constructor
    //TODO: add default ParamDict
    System(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~System();

    //Make changes to System state
    void do_paramdict_assign(ParamDict &theParams);
    void do_particle_init();
    void do_spring_init();
    void apply_pbc(); //apply periodic boundary conditions
    void zero_com(); //subtract center of mass from each particle position
    void set_obs(Observer &anObs);

    //Get (some of these could be made static)
    arma::vec get_com();
    double get_energy();
    double get_energy_between(Particle &p1, Particle &p2);
    double get_energy_cell_list();
    double get_energy_neighbor_grid();
    double get_lj_potential(double r, double sig, double eps, double rc);
    double get_wca_potential(double r, double sig, double eps);
    double get_harmonic_potential(double r, double KK, double ll0);
    double get_fene_potential(double r, double KK, double ll0, double dr);
    std::vector<arma::vec> get_forces();
    std::vector<arma::vec> get_forces_cell_list();
    arma::vec get_force(Particle &p1);
    arma::vec get_force_from(Particle &p1, Particle &p2);
    arma::vec get_force_neighbor_grid(Particle &p);
    arma::vec get_lj_force(double r, arma::vec rvec, double sig, double eps);
    arma::vec get_wca_force(double r, arma::vec rvec, double sig, double eps);
    arma::vec get_harmonic_force(double r, arma::vec rvec, double K, double l0);
    arma::vec get_fene_force(double r, arma::vec rvec, double K, double l0, double dr);
    double minimize_energy(double tol=1e-6);
    arma::vec get_disp_vec(Particle &p1, Particle &p2);
    double get_dist(Particle &p1, Particle &p2);
    Observer get_obs();

    //Methods for cell list and neighbor grid
    int needs_update_cell_list();
    void create_cell_list();
    void fill_cellneigh();
    void update_neighborgrid();
};

#endif