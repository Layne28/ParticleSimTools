#include "ActiveSolver.hpp"

ActiveSolver::ActiveSolver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg) : Solver(theParams, the_rg)
{
    if(theParams.is_key("va")) va = std::stod(theParams.get_value("va"));
    //Modify parameters to make particles and active noise generator consistent 
    theParams.add_entry("D", std::to_string(va*va)); //Active noise generator will generate values with variance va^2
    double dee_x = theSys.edges[0]/int(std::stoi(theParams.get_value("nx")));
    //Set precision of dx
    std::ostringstream out;
    out.precision(15);
    out << std::fixed << dee_x;
    theParams.add_entry("dx", out.str());
    anGen = new Generator(theParams, the_rg);

    std::cout << anGen->ny << std::endl;

    //Check that active noise gen. params are consistent with particle system
    if(fabs(anGen->Lx-theSys.edges[0])>1e-10) {
        std::cout << "ERROR! System and Active Noise Generator do not have the same size (" << std::fixed << theSys.edges[0] << " vs " << anGen->Lx << ")." << std::endl;
        exit(1);
    }
}

ActiveSolver::~ActiveSolver()
{
    delete anGen;
}

void ActiveSolver::update(System &theSys)
{
    //Get conservative forces
    std::vector<arma::vec> potential_forces = theSys.get_forces();

    //Get active noise force on each particle
    std::vector<arma::vec> active_forces;
    active_forces = get_active_noise_forces(theSys, *anGen);
    
    for (int i=0; i<theSys.N; i++) {
        theSys.particles[i].old_pos = theSys.particles[i].pos;
        arma::vec pos = theSys.particles[i].get_pos();
        for (int k=0; k<theSys.dim; k++) {
            if (fabs(potential_forces[i][k])>10000.0) std::cout << "Warning: large lj force: " << potential_forces[i][k] << std::endl;
            //Euler step
            double incr = potential_forces[i](k)/gamma*dt 
                        + active_forces[i](k)*dt;
            if (theSys.kT>0) incr += sqrt(2*theSys.kT/gamma)*gsl_ran_gaussian(rg, sqrt(dt)); 
            pos(k) += incr; 
                   
            theSys.particles[i].pos[k] = pos(k);
            theSys.particles[i].vel[k] = incr/dt;
        }
    }
    anGen->step(dt); //advance active noise in time
    theSys.apply_pbc();
    theSys.time++;
    if (theSys.do_neighbor_grid) theSys.update_neighborgrid();
}

void ActiveSolver::update_adaptive(System &theSys, double deet, int level){

    //Get conservative forces
    std::vector<arma::vec> potential_forces = theSys.get_forces();

    //Get active noise force on each particle
    std::vector<arma::vec> active_forces;
    active_forces = get_active_noise_forces(theSys, *anGen);

    //Get thermal force on each particle
    std::vector<arma::vec> thermal_forces;
    thermal_forces = get_thermal_forces(theSys, deet);
    
    //Compute particle motion due to forces
    std::vector<arma::vec> incr(theSys.N);
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        incr[i] = v;
    }
    for(int i=0; i<theSys.N; i++){
        incr[i] = potential_forces[i]/gamma*deet 
                    + active_forces[i]*deet
                    + thermal_forces[i];   
    }

    //Check whether the new position will result in a really large force
    for(int i=0; i<theSys.N; i++){
        theSys.particles[i].old_pos = theSys.particles[i].pos;
        theSys.particles[i].pos += incr[i];
    }
    theSys.apply_pbc();
    std::vector<arma::vec> new_forces = theSys.get_forces();
    double max_force = 0;
    for(int i=0; i<theSys.N; i++){
        for(int k=0; k<theSys.dim; k++){
            if(fabs(new_forces[i][k])>max_force) max_force = new_forces[i][k];
        }
    }
    //only decrease time step if force is above threshold
    //and timestep is not already tiny
    if(max_force > force_thresh && deet>1e-6){
        std::cout << "Force too high. Decreasing time step by a factor of 4 (now =" << deet/4 << ")." << std::endl;
        //revert to old position
        for(int i=0; i<theSys.N; i++){
            theSys.particles[i].pos = theSys.particles[i].old_pos;
        }
        for(int k=0; k<4; k++){
            update_adaptive(theSys, deet/4, level+1);
        }
    }
    else{
        for(int i=0; i<theSys.N; i++){
            theSys.particles[i].vel = incr[i]/deet;
        }
        anGen->step(deet);
        theSys.apply_pbc();
        if (theSys.do_neighbor_grid) theSys.update_neighborgrid();    
    }

    if(level==0) theSys.time++;

}

std::vector<arma::vec> ActiveSolver::get_active_noise_forces(System &theSys, Generator &gen)
{
    std::vector<arma::vec> active_forces(theSys.N);

    //Initialize to zero
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        active_forces[i] = v;
    } 

    //if active velocity is zero, don't bother doing calculation
    if (va<1e-10) {
        //std::cout << "Active speed is zero." << std::endl;
        return active_forces;
    }

    //Compute real-space noise
    arma::field<arma::cx_vec> xi = gen.get_xi_r(1); //do_fft=1

    if (theSys.dim==3) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = int((pos(0)+0.5*theSys.edges[0])/gen.dx);
            int yind = int((pos(1)+0.5*theSys.edges[1])/gen.dx);
            int zind = int((pos(2)+0.5*theSys.edges[2])/gen.dx);
            for(int j=0; j<3; j++) active_forces[i](j) = xi(xind,yind,zind)(j).real(); 
        }
    }
    else if (theSys.dim==2) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = int((pos(0)+0.5*theSys.edges[0])/gen.dx);
            int yind = int((pos(1)+0.5*theSys.edges[1])/gen.dx);
            for(int j=0; j<2; j++) active_forces[i](j) = xi(xind,yind)(j).real(); 
            if (fabs(active_forces[i][0])>100.0) std::cout << "Warning: large active force: " << active_forces[i][0] << std::endl;
            if (fabs(active_forces[i][1])>100.0) std::cout << "Warning: large active force: " << active_forces[i][1] << std::endl;
        }
    }
    else if (theSys.dim==1) {
        //Assign active force to each particle based on location in noise grid
        for(int i=0; i<theSys.N; i++) {
            arma::vec pos = theSys.particles[i].get_pos();
            int xind = int((pos(0)+0.5*theSys.edges[0])/gen.dx);
            active_forces[i](0) = xi(xind)(0).real(); 
        }
    }

    return active_forces;
}

std::vector<arma::vec> ActiveSolver::get_thermal_forces(System &theSys, double deet){
    
    std::vector<arma::vec> thermal_forces(theSys.N);

    //Initialize to zero
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        thermal_forces[i] = v;
    }

    //if temperature is zero, don't bother doing calculation
    if (theSys.kT<1e-10) {
        //std::cout << "kT is zero." << std::endl;
        return thermal_forces;
    }

    for (int i=0; i<theSys.N; i++) {
        for (int k=0; k<theSys.dim; k++) {
            thermal_forces[i][k] = sqrt(2*theSys.kT/gamma)*gsl_ran_gaussian(rg, sqrt(deet)); 
        }
    }

    return thermal_forces;
}