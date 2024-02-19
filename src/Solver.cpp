#include "Solver.hpp"

Solver::Solver(System &theSys, ParamDict &theParams, gsl_rng *&the_rg)
{
    dt = 0.005;
    gamma = 1.0;

    if(theParams.is_key("dt")) dt = std::stod(theParams.get_value("dt"));
    if(theParams.is_key("gamma")) gamma = std::stod(theParams.get_value("gamma"));
    if(theParams.is_key("do_adaptive_timestep")) do_adaptive_timestep = std::stoi(theParams.get_value("do_adaptive_timestep")); 
    if(theParams.is_key("do_subtract_com")) do_subtract_com = std::stoi(theParams.get_value("do_subtract_com")); 

    //***Active noise***
    if(theParams.is_key("va")) va = std::stod(theParams.get_value("va"));
    if(theParams.is_key("nx")) nx = std::stoi(theParams.get_value("nx")); 
    if(theParams.is_key("ny")) ny = std::stoi(theParams.get_value("ny")); 
    if(theParams.is_key("nz")) nz = std::stoi(theParams.get_value("nz")); 
    if(theParams.is_key("interpolation")) interpolation = theParams.get_value("interpolation"); 
    //Active noise generator will generate values with variance va^2
    theParams.add_entry("D", std::to_string(va*va)); 

    //Set grid spacing
    double dee_x = theSys.edges[0]/nx;
    //Set precision of dx
    std::ostringstream out_x;
    out_x.precision(15);
    out_x << std::fixed << dee_x;
    theParams.add_entry("dx", out_x.str());

    if(theSys.dim==2 || theSys.dim==3){
        double dee_y = theSys.edges[1]/ny;
        //Set precision of dy
        std::ostringstream out_y;
        out_y.precision(15);
        out_y << std::fixed << dee_y;
        theParams.add_entry("dy", out_y.str());
    }

    if(theSys.dim==3){
        double dee_z = theSys.edges[2]/nz;
        //Set precision of dz
        std::ostringstream out_z;
        out_z.precision(15);
        out_z << std::fixed << dee_z;
        theParams.add_entry("dz", out_z.str());
    }
    
    anGen = std::make_unique<Generator>(theParams, the_rg);
    //Check that active noise gen. params are consistent with particle system
    if(fabs(anGen->Lx-theSys.edges[0])>1e-10) {
        std::cout << "ERROR! System and Active Noise Generator do not have the same size (" << std::fixed << theSys.edges[0] << " vs " << anGen->Lx << ")." << std::endl;
        exit(1);
    }

    //***Set RNG***
    rg = the_rg;
}

Solver::~Solver() 
{
}

void Solver::update(System &theSys, double deet, int level)
{
    //Get conservative forces
    std::vector<arma::vec> potential_forces = theSys.get_forces();

    //Get AOUP forces
    std::vector<arma::vec> aoup_forces = get_aoup_forces(theSys);

    //Get active noise forces
    std::vector<arma::vec> active_forces = get_active_noise_forces(theSys, *anGen);

    //Get thermal forces
    std::vector<arma::vec> thermal_forces = get_thermal_forces(theSys, deet);

    //Compute particle motion due to forces
    std::vector<arma::vec> incr(theSys.N);
    for(int i=0; i<theSys.N; i++){
        arma::vec v(theSys.dim,arma::fill::zeros);
        incr[i] = v;
    }
    for(int i=0; i<theSys.N; i++){
        incr[i] = potential_forces[i]/gamma*deet
                  + aoup_forces[i]/gamma*deet
                  + active_forces[i]*deet //TODO: check whether there should be a factor of 1/gamma here
                  + thermal_forces[i]; //deet is included in thermal forces
    }

    //Update positions accordingly
    for(int i=0; i<theSys.N; i++){
        theSys.particles[i].old_pos = theSys.particles[i].pos;
        theSys.particles[i].pos += incr[i];
    }
    theSys.apply_pbc();

    if(do_adaptive_timestep){
        //Check whether the new position will result in a really large force
        //TODO: How to modify this for bond breaking?
        std::vector<arma::vec> new_forces = theSys.get_forces();
        double max_force = 0;
        for(int i=0; i<theSys.N; i++){
            for(int k=0; k<theSys.dim; k++){
                if(fabs(new_forces[i][k])>max_force) max_force = new_forces[i][k];
            }
        }
        //Only decrease time step if force is above threshold
        //and timestep is not already tiny
        if(max_force > force_thresh && deet>1e-10){
            if(level==0) std::cout << "Force too high. Decreasing time step by a factor of 4 (now =" << deet/4 << ")." << std::endl;
            //revert to old position
            for(int i=0; i<theSys.N; i++){
                //TODO: WARNING: to make this work with NeighborGrid, need to reset old_pos as well.
                theSys.particles[i].pos = theSys.particles[i].old_pos;
            }
            for(int k=0; k<4; k++){
                update(theSys, deet/4, level+1);
            }
        }
        else{
            for(int i=0; i<theSys.N; i++){
                theSys.particles[i].vel = incr[i]/deet;
                //Advance OU process in time
                if(theSys.particles[i].is_aoup){
                    update_self_prop_vel(theSys, i, deet);
                }
            }
            if(va>0.0) anGen->step(deet); //advance active noise in time
            theSys.apply_pbc();
        }
    }
    else{
        for(int i=0; i<theSys.N; i++){
            theSys.particles[i].vel = incr[i]/deet;
            //Advance OU process in time
            if(theSys.particles[i].is_aoup){
                update_self_prop_vel(theSys, i, deet);
            }
        }
        if(va>0.0) anGen->step(deet); //advance active noise in time
        theSys.apply_pbc();
    }

    //If applicable, break and make bonds
    if(theSys.can_bonds_break==1 && level==0){
        update_bonds(theSys, deet);
    }

    if(level==0) theSys.time++;
}

void Solver::update_bonds(System &theSys, double deet){

    double k0 = theSys.k0_bond;

    //Keep track of which bonds have been broken
    //so that we don't attempt to remake them this timestep
    std::vector<std::pair<int, int>> broken_bonds;

    //Attempt to break bonds
    double P_detach = 1.0 - exp(-k0*deet);
    for(int i=0; i<theSys.N; i++){
        int nsprings = theSys.particles[i].get_num_springs();
        std::vector<int> to_remove;
        for(int j=0; j<nsprings; j++){
            //To prevent double-counting, only attempt a bond break if 
            //the second particle index is larger than the first
            Particle *p2 = theSys.particles[i].springs[j].node2;
            if(theSys.particles[i].is_equal(*p2)) p2 = theSys.particles[i].springs[j].node1;
            if(p2->get_id()>theSys.particles[i].get_id()){
                double xsi = gsl_rng_uniform(rg);
                if (xsi<P_detach){
                    to_remove.push_back(j);
                    //Record broken bond
                    int pid1 = theSys.particles[i].get_id();
                    int pid2 = p2->get_id();
                    std::pair<int, int> bond = std::make_pair(pid1, pid2);
                    broken_bonds.push_back(bond);
                }
            }
        }
        //sort lowest to highest spring index
        std::sort(to_remove.begin(), to_remove.end());
        //Remove bonds starting from the highest index
        //so that we don't have to worry about re-indexing
        while(!to_remove.empty()){
            Particle *p2 = theSys.particles[i].springs[to_remove.back()].node2;
            if(theSys.particles[i].is_equal(*p2)) p2 = theSys.particles[i].springs[to_remove.back()].node1;
            Spring::remove_spring(theSys.particles[i], *p2);
            //std::cout << "Broke a bond!" << std::endl;
            to_remove.pop_back();
        }
    }

    //Attempt to make bonds
    //Don't try making bonds b/t particles further away than 
    //max bond extension
    double max_dist = theSys.l0 + theSys.drmax;

    if(theSys.do_cell_list){
        //update cell list if particles have moved out of cells
        if(theSys.needs_update_cell_list()==1){
            theSys.create_cell_list();
        }
        //Loop over cells (including self)
        for (int index1 = theSys.N-1; index1 >= 0; index1--) {
            int icell = theSys.cellndx[index1];
            int index2 = theSys.head[icell];
            for (int nc = 0; nc < theSys.cellneigh[icell][0]; nc++) {
                int jcell = theSys.cellneigh[icell][nc+1];
                index2 = theSys.head[jcell];
                while (index2 != -1) {
                    if(index1>index2){
                        int pid1 = theSys.particles[index1].get_id();
                        int pid2 = theSys.particles[index2].get_id();
                        std::pair<int,int> bond = {pid1, pid2};
                        //Order small to large, same as broken bonds list
                        if(pid2<pid1){
                            bond = {pid2, pid1};
                        }
                        //Check if particles are close but not bonded
                        if(theSys.get_dist(theSys.particles[index1], theSys.particles[index2])<max_dist 
                            && !(theSys.particles[index1].has_connection(theSys.particles[index2]))){
                            //Check that this is not a just-broken bond
                            if(std::find(broken_bonds.begin(), broken_bonds.end(), bond)==broken_bonds.end()){
                                double bond_energy = theSys.get_bonded_energy(theSys.particles[index1], theSys.particles[index2]);
                                double k_on = k0*std::exp(-bond_energy/theSys.kT);
                                double P_attach = 1 - std::exp(-k_on*deet);
                                double xsi = gsl_rng_uniform(rg);
                                if (xsi<P_attach){
                                    Spring::add_spring(theSys.particles[index1], theSys.particles[index2]);
                                }
                            }
                        }
                    }
                    index2 = theSys.list[index2];
                }
            }
        }
    }

    //Do O(N^2) loop
    else{
        for(int i=0; i<theSys.N-1; i++){
            for(int j=i+1; j<theSys.N; j++){
                int pid1 = theSys.particles[i].get_id();
                int pid2 = theSys.particles[j].get_id();
                if(pid1==pid2) std::cout << "ERROR: DIFFERENT PARTICLE IDS ARE THE SAME" << std::endl;
                if(i==j) std::cout << "ERROR" << std::endl;
                std::pair<int,int> bond = {pid1, pid2};
                //Order small to large, same as broken bonds list
                if(pid2<pid1){
                    bond = {pid2, pid1};
                }
                //Check if particles are close but not bonded
                if(theSys.get_dist(theSys.particles[i], theSys.particles[j])<max_dist 
                    && !(theSys.particles[i].has_connection(theSys.particles[j]))){
                    //Check that this is not a just-broken bond
                    if(std::find(broken_bonds.begin(), broken_bonds.end(), bond)==broken_bonds.end()){
                        double bond_energy = theSys.get_bonded_energy(theSys.particles[i], theSys.particles[j]);
                        double k_on = k0*std::exp(-bond_energy/theSys.kT);
                        double P_attach = 1 - std::exp(-k_on*deet);
                        double xsi = gsl_rng_uniform(rg);
                        if (xsi<P_attach){
                            Spring::add_spring(theSys.particles[i], theSys.particles[j]);
                        }
                    }
                }
            }
        }
    }
    
}

std::vector<arma::vec> Solver::get_thermal_forces(System &theSys, double deet){
    
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

    if(do_subtract_com==1) thermal_forces = get_subtracted_com_forces(theSys, thermal_forces);

    return thermal_forces;
}

std::vector<arma::vec> Solver::get_aoup_forces(System &theSys){
    
    std::vector<arma::vec> aoup_forces(theSys.N);

    //Initialize to zero
    for (int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        aoup_forces[i] = v;
    }

    for (int i=0; i<theSys.N; i++) {
        //If D0 is zero, don't bother doing calculation
        if (theSys.particles[i].aoup_D0<1e-10){
            for(int k=0; k<theSys.dim; k++){
                aoup_forces[i][k] = 0.0;
            }
        }
        else{
            for (int k=0; k<theSys.dim; k++) {
                //TODO: Need to multiply by appropriate factors of D0, tau
                aoup_forces[i][k] = theSys.particles[i].self_prop_vel[k]; 
            }
        }
    }

    if(do_subtract_com==1) aoup_forces = get_subtracted_com_forces(theSys, aoup_forces);

    return aoup_forces;
}

std::vector<arma::vec> Solver::get_active_noise_forces(System &theSys, Generator &gen)
{
    std::vector<arma::vec> active_forces(theSys.N);

    //Initialize to zero
    for(int i=0; i<theSys.N; i++) {
        arma::vec v(theSys.dim,arma::fill::zeros);
        active_forces[i] = v;
    }

    //if active velocity is zero, don't bother doing calculation
    if (va<1e-10) {
        return active_forces;
    }

    //Compute real-space noise
    arma::field<arma::cx_vec> xi = gen.get_xi_r(1);
    
    //Assign active force to each particle
    //based on location in noise grid
    //TODO: write tests for this
    for(int i=0; i<theSys.N; i++)
    {
        arma::vec pos = theSys.particles[i].pos;
        std::vector<int> indices(theSys.dim, 0);
        for(int k=0; k<theSys.dim; k++){
            indices[k] = int((pos(k)+0.5*theSys.edges[k])/gen.spacing[k]);
        }
        for(int k=0; k<theSys.dim; k++){
            if(theSys.dim==1){
                if (interpolation=="none"){
                    active_forces[i](k) = xi(indices[0])(k).real();
                }
                else if (interpolation=="linear"){
                    //Treat noise as being at midpoint xmi of cell i.
                    //If particle has position xj, xmi<=xj<xmi+dx/2,
                    //then compute noise at x as: 
                    //xi(x) = l*xi(xm(i+1)) + (1-l)*xi(xmi), l=(x-xmi)/dx
                    //If particle has position xj, xmi-dx/2<=xj<xmi,
                    //then compute noise at x as: 
                    //xi(x) = l*xi(xm(i+1)) + (1-l)*xi(xmi), l=(x-xmi)/dx
                    double scaled_pos = (pos(k)+0.5*theSys.edges[k])/gen.spacing[k];
                    double whole_pos = floor(scaled_pos)
                    double decimal_pos = scaled_pos - whole_pos;
                    if(decimal_pos)>=0.5{
                        int ind = int((pos(k)+0.5*theSys.edges[k])/gen.spacing[k]);
                        double xi1 = xi(ind)(k).real();
                        double xi2 = xi((ind+1)%nx)(k).real();
                        double ell = scaled_pos - (whole_pos + 0.5);
                        double interp_xi = ell*xi2 + (1-ell)*xi1;
                    }
                    else{
                        int ind = int((pos(k)+0.5*theSys.edges[k])/gen.spacing[k]);
                        double xi1 = xi(ind)(k).real();
                        double xi0 = xi((ind-1+nx)%nx)(k).real();
                        double ell = scaled_pos - (whole_pos - 0.5);
                        double interp_xi = ell*xi1 + (1-ell)*xi0;
                    }
                }
                else{
                    std::cout << "ERROR: interpolation method not recognized." << std::endl;
                    exit(-1);
                }
            }
            else if(theSys.dim==2){
                if (interpolation=="none"){
                    active_forces[i](k) = xi(indices[0],indices[1])(k).real();
                }
                else{
                    std::cout << "ERROR: interpolation method not recognized." << std::endl;
                    exit(-1);
                }
            }
            else{
                if (interpolation=="none"){
                    active_forces[i](k) = xi(indices[0],indices[1],indices[2])(k).real();
                }
                else{
                    std::cout << "ERROR: interpolation method not recognized." << std::endl;
                    exit(-1);
                }
            } 
        } 
    }

    if(do_subtract_com==1) active_forces = get_subtracted_com_forces(theSys, active_forces);

    return active_forces;
}

std::vector<arma::vec> Solver::get_subtracted_com_forces(System &theSys, std::vector<arma::vec> forces){

    arma::vec com_vel(theSys.dim, arma::fill::zeros);
    for(int i=0; i<theSys.N; i++) {
        com_vel += forces[i];
    }
    for(int i=0; i<theSys.N; i++) {
        forces[i] -= com_vel/theSys.N;
    }

    return forces;
}

void Solver::update_self_prop_vel(System &theSys, int index, double deet){

    double D0 = theSys.particles[index].aoup_D0;
    double tau = theSys.particles[index].aoup_tau;
    for(int k=0; k<theSys.dim; k++){
        theSys.particles[index].self_prop_vel[k] += (-theSys.particles[index].self_prop_vel[k]/tau)*deet + sqrt(2*D0)/tau*gsl_ran_gaussian(rg, sqrt(deet));
        //std::cout << theSys.particles[index].self_prop_vel(k) << std::endl;
    }
}