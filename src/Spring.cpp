#include <iostream>
#include <cmath>
#include "Spring.hpp"
#include "Particle.hpp" 

Spring::Spring(Particle& n1, Particle& n2, double K_1, double l0_1)
{
    K = K_1;
    l0 = l0_1;
    node1 = &n1;
    node2 = &n2;
}

Spring::~Spring()
{
    //TODO: Delete pointers?
}

void Spring::set_stiffness(double K_1)
{
    this->K = K_1;
}

void Spring::set_rest_length(double l0_1)
{
    this->l0 = l0_1;
}

double Spring::get_stiffness()
{
    return this->K;
}

double Spring::get_rest_length()
{
    return this->l0;
}

double Spring::get_length()
{
    arma::vec pos1 = (*(this->node1)).get_pos();
    arma::vec pos2 = (*(this->node2)).get_pos();
    double len = 0;
    for (int i=0; i<3; i++)
    {
        len += (pos1[i]-pos2[i])*(pos1[i]-pos2[i]);
    }
    len = sqrt(len);

    return len;
}

double Spring::get_strain()
{
    double ell = this->get_length();
    double ell0 = this->get_rest_length();

    return (ell-ell0)/ell0;
}

bool Spring::is_spring(Particle &n1, Particle &n2)
{
    return (n1.has_connection(n2) && n2.has_connection(n1));
}

Spring Spring::add_spring(Particle &n1, Particle &n2, double K_1, double l0_1)
{
    
    if (Spring::is_spring(n1, n2))
    {
        throw std::runtime_error("Cannot create Spring -- Particles are already connected!");
    }
    
    Spring s(n1, n2, K_1, l0_1);
    n1.springs.push_back(s);
    n2.springs.push_back(s);

    return s;
}

void Spring::remove_spring(Particle &n1, Particle &n2)
{
    if (!Spring::is_spring(n1, n2))
    {
        throw std::runtime_error("Cannot remove Spring -- it doesn't exist!");
    }
    for (int i=0; i<n1.get_num_springs(); i++)
    {
        if (n1.springs[i].node1->is_equal(n2) || n1.springs[i].node2->is_equal(n2))
        {
            n1.springs.erase(n1.springs.begin()+i);
            break;
        }
    }
    for (int i=0; i<n2.get_num_springs(); i++)
    {
        if (n2.springs[i].node1->is_equal(n1) || n2.springs[i].node2->is_equal(n1))
        {
            n2.springs.erase(n2.springs.begin()+i);
            return;
        }
    }
    throw std::runtime_error("Error -- could not find Spring to remove!");
}

Spring Spring::get_spring(Particle &n1, Particle &n2)
{
    if (!Spring::is_spring(n1, n2))
    {
        throw std::runtime_error("Cannot get Spring -- it doesn't exist!");
    }

    for (int i=0; i<n1.get_num_springs(); i++)
    {
        if (n1.springs[i].node1->is_equal(n2) || n1.springs[i].node2->is_equal(n2))
        {
            return n1.springs[i];
        }
    }
    for (int i=0; i<n2.get_num_springs(); i++)
    {
        if (n2.springs[i].node1->is_equal(n1) || n2.springs[i].node2->is_equal(n1))
        {
            return n2.springs[i];
        }
    }
    throw std::runtime_error("Error -- could not find Spring!");
}

std::ostream& operator<<(std::ostream& os, Spring& s)
{

    return os << "Spring with Stiffness: " << s.get_stiffness() << "; Rest Length: " << s.get_rest_length() << "; Particles:\n" << "    " << *s.node1 << "\n" << "    " << *s.node2 << std::endl;
}