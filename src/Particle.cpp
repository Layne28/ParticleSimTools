#include <iostream>
#include "Particle.hpp"

//TODO: write unit tests for object id/counter
int Particle::counter = 0;

Particle::Particle(int dim, int is_a_node, int is_an_aoup, double my_aoup_D0, double my_aoup_tau) {
    d = dim;
    pos.zeros(dim);
    old_pos.zeros(dim);
    vel.zeros(dim);
    active_force.zeros(dim);
    conservative_force.zeros(dim);
    id = Particle::counter;
    Particle::counter++;

    is_node = is_a_node;
    is_aoup = is_an_aoup;
    if(is_aoup==1){
        self_prop_vel.zeros(dim);
        aoup_D0 = my_aoup_D0;
        aoup_tau = my_aoup_tau;
    }
}

Particle::Particle(int dim, int is_a_node) {
    d = dim;
    pos.zeros(dim);
    old_pos.zeros(dim);
    vel.zeros(dim);
    active_force.zeros(dim);
    conservative_force.zeros(dim);
    id = Particle::counter;
    Particle::counter++;

    is_node = is_a_node;
    is_aoup = 0;
    if(is_aoup==1){
        self_prop_vel.zeros(dim);
        aoup_D0 = 0.0;
        aoup_tau = 1.0;
    }
}

Particle::Particle(int dim) {
    d = dim;
    pos.zeros(dim);
    old_pos.zeros(dim);
    vel.zeros(dim);
    active_force.zeros(dim);
    conservative_force.zeros(dim);
    id = Particle::counter;
    Particle::counter++;

    is_node = 0;
    is_aoup = 0;
    if(is_aoup==1){
        self_prop_vel.zeros(dim);
        aoup_D0 = 0.0;
        aoup_tau = 1.0;
    }
}

Particle::~Particle() {}

int Particle::get_id() {
    return id;
}

int Particle::get_num_springs()
{
    if(!this->is_node) return 0;
    else return this->springs.size();
}

arma::vec Particle::get_pos() {
    return pos;
}

void Particle::set_pos(arma::vec mypos) {
    if(mypos.n_elem!=this->d){
        std::cout << "Error: attempting to assign position with "
        " wrong dimensonal vector." << std::endl;
        exit(-1);
    }
    for(int k=0; k<this->d; k++){
        this->pos[k] = mypos[k];
    }
}

bool Particle::is_equal(Particle &p){
    if (this->id==p.id) {
        return true;
    }
    else {
        return false;
    }
}

bool Particle::has_connection(Particle &n)
{
    if (!n.is_node) return false;
    else{
        int nspring = this->get_num_springs();
        for (int i=0; i<nspring; i++)
        {
            if (this->springs[i].node1->is_equal(n) || this->springs[i].node2->is_equal(n))
            {
                return true;
            }
        }
        return false;
    }
}

std::ostream& operator<<(std::ostream& os, Particle& p) {
    arma::vec pos = p.get_pos();
    if (p.d==3){
        return os << "Particle with id: " << p.get_id() << " and Pos: (" << pos(0) << ", " << pos(1) << ", " << pos(2) << ")";
    } else if(p.d==2){
       return os << "Particle with id: " << p.get_id() << " and Pos: (" << pos(0) << ", " << pos(1) << ")"; 
    }
    else{
        std::cout << "Error: num of dims not supported." << std::endl;
        exit(-1);
    }
}

double* Particle::get_position() {
    return pos.memptr();
}

double* Particle::get_old_position() {
    return old_pos.memptr();
}
