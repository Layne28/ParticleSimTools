#include <iostream>
#include "Particle.hpp"

//TODO: write unit tests for object id/counter
int Particle::counter = 0;

Particle::Particle(int dim) {
    d = dim;
    pos.zeros(dim);
    old_pos.zeros(dim);
    vel.zeros(dim);
    id = Particle::counter;
    Particle::counter++;
}

Particle::~Particle() {}


int Particle::get_id() {
    return id;
}

arma::vec Particle::get_pos() {
    return pos;
}

bool Particle::is_equal(Particle &p){
    if (this->id==p.id) {
        return true;
    }
    else {
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
