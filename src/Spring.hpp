//A Spring object has a rest length, a stiffness, and pointers to the two Particles it connects.

#ifndef SPRING_HPP
#define SPRING_HPP

class Particle;

class Spring
{
private:
    friend class Particle;

public:
    double K; //stiffness
    double l0; //rest length
    Particle* node1;
    Particle* node2;

    //constructor (make this private in the future?)
    Spring(Particle& n1, Particle& n2, double K_1=0.0, double l0_1=0.0);

    //destructor
    ~Spring();

    void set_stiffness(double K_1);
    void set_rest_length(double l0_1);
    double get_stiffness();
    double get_rest_length();

    //Move these two to System
    double get_length(); //need to make this respect pbc
    double get_strain(); //need to make this respect pbc

    static bool is_spring(Particle &n1, Particle &n2);
    static Spring add_spring(Particle &n1, Particle &n2, double K_1=0.0, double l0_1=0.0);
    static void remove_spring(Particle &n1, Particle &n2);
    static Spring get_spring(Particle &n1, Particle &n2);

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, Spring& s);
};

#endif