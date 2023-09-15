#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <array>
#include <vector>
#include <armadillo>
#include "Spring.hpp"

class Particle
{
private:
    static int counter; //Warning -- not thread-safe!!
    int id;

    friend class Spring;
    
public:
    int d;
    arma::vec pos;
    arma::vec old_pos;
    arma::vec vel;

    //If the particle is part of a network
    int is_node;
    std::vector<Spring> springs;

    //If the particle is an AOUP
    int is_aoup;
    double aoup_D0;
    double aoup_tau;
    arma::vec self_prop_vel;

    //constructor
    Particle(int dim, int is_an_aoup=0, double my_aoup_D0=0.0, double my_aoup_tau=0.0);

    //destructor
    ~Particle();

    //Various functions
    int get_id();
    arma::vec get_pos();
    bool is_equal(Particle &p);

    bool has_connection(Particle &p);
    int get_num_springs();
    void move();

    arma::vec get_self_prop_vel();

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, Particle& p);

    double *get_position();
    double *get_old_position();
};

#endif