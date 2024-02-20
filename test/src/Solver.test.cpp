#include "catch.hpp"
#include "../../src/Solver.hpp"
#include "../../src/Particle.hpp"
#include "../../src/Spring.hpp"
#include "../../src/System.hpp"
#include "../../src/CustomRandom.hpp"

TEST_CASE("Interpolation"){
    ParamDict defaultParams;
    defaultParams.add_entry("particle_protocol", "zeros");
    defaultParams.add_entry("dim", "1");
    defaultParams.add_entry("phi", "0.0");
    defaultParams.add_entry("Lx", "3.0");
    defaultParams.add_entry("nx", "6");
    defaultParams.add_entry("is_p_x", "1");
    defaultParams.add_entry("do_cell_list", "0");
    defaultParams.add_entry("interpolation", "linear");

    gsl_rng *myGen = CustomRandom::init_rng(1);
    System theSys(defaultParams, myGen);
    Solver solver(theSys, defaultParams, myGen);

    SECTION("Check 1d interpolation function"){
        double x = 0.7;
        double f1 = 3.0;
        double f2 = 4.0;
        double interp = solver.interpolate_1d(x, f1, f2);
        REQUIRE(std::fabs(interp-3.7)<1e-10);
    }

    SECTION("Check 1d interpolation of noise"){
        theSys.particles[0].pos[0] = 0.0-theSys.edges[0]/2.0;
        std::vector<arma::vec> active_noise = solver.get_active_noise_forces(theSys, *solver.anGen);
        arma::field<arma::cx_vec> xi = (*solver.anGen).get_xi_r(1);

        //Do interpolation by hand to test
        double thenoise = xi(solver.nx-1)(0).real()*0.5 + xi(0)(0).real()*0.5;
        REQUIRE(std::fabs(active_noise[0](0)-thenoise)<1e-10);

        //Try another position
        theSys.particles[0].pos[0] = 0.1-theSys.edges[0]/2.0;
        active_noise = solver.get_active_noise_forces(theSys, *solver.anGen);
        thenoise = xi(solver.nx-1)(0).real()*0.3 + xi(0)(0).real()*0.7;
        REQUIRE(std::fabs(active_noise[0](0)-thenoise)<1e-10);

        //One more
        theSys.particles[0].pos[0] = 2.35-theSys.edges[0]/2.0;
        active_noise = solver.get_active_noise_forces(theSys, *solver.anGen);
        thenoise = xi(solver.nx-2)(0).real()*0.8 + xi(solver.nx-1)(0).real()*0.2;
        REQUIRE(std::fabs(active_noise[0](0)-thenoise)<1e-10);
    }

    /*
    SECTION("Check greater than midpoint"){
        double x = 1.6;
        double f0 = 2.0;
        double f1 = 3.0;
        double f2 = 4.0;
        double interp = solver.interpolate_1d(x, f0, f1, f2);
        REQUIRE(std::fabs(interp-3.1)<1e-10);
    }
    */

}

TEST_CASE("Solver")
{
    ParamDict defaultParams;
    gsl_rng *myGen = CustomRandom::init_rng(1);
    defaultParams.add_entry("k0_bond", "1.0");
    defaultParams.add_entry("kT", "1.0");
    defaultParams.add_entry("particle_protocol", "uniform_lattice");
    defaultParams.add_entry("dim", "1");
    defaultParams.add_entry("is_p_x", "0");
    defaultParams.add_entry("Lx", "2.0");
    defaultParams.add_entry("Nx", "2");
    defaultParams.add_entry("a", "1");
    defaultParams.add_entry("dt", "0.01");
    defaultParams.add_entry("is_network", "1");
    defaultParams.add_entry("do_cell_list", "0");

    SECTION("Check single bond detachment")
    {
        //Generate bond lifetimes for different realizations
        int nsim=10000;
        double lifetimes[nsim] = {0};
        for(int s=0; s<nsim; s++){
            System theSys(defaultParams, myGen);
            Solver solver(theSys, defaultParams, myGen);

            //Only one check should be necessary here
            if(s==0){
                REQUIRE(theSys.get_dist(theSys.particles[0],theSys.particles[1])==1.0);
                REQUIRE(theSys.particles[0].get_num_springs()==theSys.particles[1].get_num_springs());
            }
            int is_broken = 0;
            while(is_broken==0){
                solver.update_bonds(theSys, theSys.dt);
                if(s==0) REQUIRE(theSys.particles[0].get_num_springs()==theSys.particles[1].get_num_springs());
                if(theSys.particles[0].get_num_springs()==0){
                    is_broken = 1;
                    lifetimes[s] = theSys.dt*theSys.time;
                }
                theSys.time++;
            }
        }

        //Create histogram of lifetimes
        double binwidth=0.1;
        int nbins = std::round(*std::max_element(lifetimes,lifetimes+nsim)/binwidth)+1;
        double histogram[nbins] = {0};
        for(int s=0; s<nsim; s++){
            int index = std::floor(lifetimes[s]/binwidth);
            histogram[index]++;
        }
        //Normalize histogram
        for(int i=0; i<nbins; i++) histogram[i] /= nsim;

        //Print results
        std::cout << "lifetime\tprobability\tpredicted" << std::endl;
        for(int i=0; i<nbins/3; i++){
            std::cout << i*binwidth << "\t" << histogram[i] << "\t" << std::exp(-(i+0.5)*binwidth)*binwidth << std::endl;
            //Check for more than than 20% relative error in bins with significant probability
            //If this fails, try increasing "nsim"
            if(histogram[i]>0.3*histogram[0]){
                REQUIRE(std::fabs(histogram[i]-std::exp(-(i+0.5)*binwidth)*binwidth)/(std::exp(-(i+0.5)*binwidth)*binwidth)<0.20); 
            }
        }
    }

    SECTION("Check single bond attachment")
    {
        //Generate time to "react"/form bond for different realizations
        int nsim=20000;
        double eps_bond;
        double dt;
        double kT;
        double k_on;
        double rxn_times[nsim] = {0};
        defaultParams.add_entry("eps_bond", "-2.0");
        for(int s=0; s<nsim; s++){
            std::cout << "Sim. " << s << std::endl;
            System theSys(defaultParams, myGen);
            eps_bond = theSys.eps_bond;
            kT = theSys.kT;
            k_on = theSys.k0_bond*std::exp(-eps_bond/kT);
            if(theSys.particles[0].has_connection(theSys.particles[1])){
                Spring::remove_spring(theSys.particles[0], theSys.particles[1]);
            }
            REQUIRE(theSys.particles[0].get_num_springs()==0);
            Solver solver(theSys, defaultParams, myGen);
            dt = theSys.dt;

            //Only one check should be necessary here
            if(s==0){
                REQUIRE(theSys.get_dist(theSys.particles[0],theSys.particles[1])==1.0);
                REQUIRE(theSys.particles[0].get_num_springs()==theSys.particles[1].get_num_springs());
            }
            int is_broken = 1;
            while(is_broken==1){
                solver.update_bonds(theSys, theSys.dt);
                if(s==0) REQUIRE(theSys.particles[0].get_num_springs()==theSys.particles[1].get_num_springs());
                if(theSys.particles[0].get_num_springs()==1){
                    is_broken = 0;
                    rxn_times[s] = theSys.dt*theSys.time;
                }
                theSys.time++;
            }
        }
        //for(int s=0; s<nsim; s++) std::cout << "Reaction time: " << rxn_times[s] << std::endl;

        //Create histogram of rxn times
        double binwidth=0.05;
        int nbins = std::round(*std::max_element(rxn_times,rxn_times+nsim)/binwidth)+1;
        double histogram[nbins] = {0};
        for(int s=0; s<nsim; s++){
            int index = std::floor(rxn_times[s]/binwidth);
            histogram[index]++;
        }
        //Normalize histogram
        for(int i=0; i<nbins; i++) histogram[i] /= nsim;

        //Print results
        std::cout << "rxn time\tprobability\tpredicted" << std::endl;
        //double testsum = 0.0;
        //double testsum2 = 0.0;
        for(int i=0; i<nbins/3; i++){
            std::cout << i*binwidth << "\t" << histogram[i] << "\t" << std::exp(-k_on*(i+0.5)*binwidth)*binwidth*k_on << std::endl;
            //Check for more than than 20% relative error in bins with significant probability
            //If this fails, try increasing "nsim"
            if(histogram[i]>0.3*histogram[0]){
                REQUIRE(std::fabs(histogram[i]-std::exp(-k_on*(i+0.5)*binwidth)*binwidth*k_on)/(std::exp(-k_on*(i+0.5)*binwidth)*binwidth*k_on)<0.20); 
            }
            //testsum += std::exp(-(i)*binwidth)*binwidth;
            //testsum2 += histogram[i];
        }
        //std::cout << "test: " << testsum << std::endl;
        //std::cout << "test2: " << testsum2 << std::endl;
    }

    SECTION("Check detailed balance")
    {
        int nsteps = 200000000;
        //Note: if k_on becomes large as here, need to 
        //decrease timestep to recover Boltzmann statistics accurately.
        defaultParams.add_entry("dt", "0.001");
        defaultParams.add_entry("eps_bond", "-3.0");
        System theSys(defaultParams, myGen);
        Solver solver(theSys, defaultParams, myGen);

        double bond_energy = theSys.get_bonded_energy(theSys.particles[0], theSys.particles[1]);

        int num_bonded = 0;
        int num_unbonded = 0;

        //Run a trajectory
        for(int t=0; t<nsteps; t++){
            solver.update_bonds(theSys, theSys.dt);
            if(theSys.particles[0].get_num_springs()==1){
                num_bonded++;
            }
            else{
                num_unbonded++;
            }
            theSys.time++;
        }
        double ratio = (1.0*num_bonded)/(1.0*num_unbonded);
        std::cout << "Bonded-to-unbonded ratio:" << std::endl;
        std::cout << "Measured\tPredicted" << std::endl;
        std::cout << ratio << "\t" << std::exp(-bond_energy/theSys.kT) << std::endl;
        double prediction = std::exp(-bond_energy/theSys.kT);
        REQUIRE(std::fabs(ratio-prediction)/prediction<1e-2);
    }

    SECTION("Check detailed balance for stretched spring")
    {
        int nsteps = 200000000;
        //Note: if k_on becomes large as here, need to 
        //decrease timestep to recover Boltzmann statistics accurately.
        defaultParams.add_entry("dt", "0.001");
        defaultParams.add_entry("l0", "1.0");
        defaultParams.add_entry("drmax", "0.5");
        defaultParams.add_entry("K", "1.0");
        defaultParams.add_entry("eps_bond", "-3.0");
        defaultParams.add_entry("bonded_potential_type", "fene");
        System theSys(defaultParams, myGen);
        Solver solver(theSys, defaultParams, myGen);

        //Stretch the spring
        theSys.particles[1].pos[0] += 0.3;

        double bond_energy = theSys.get_bonded_energy(theSys.particles[0], theSys.particles[1]);

        int num_bonded = 0;
        int num_unbonded = 0;

        //Run a trajectory
        for(int t=0; t<nsteps; t++){
            solver.update_bonds(theSys, theSys.dt);
            if(theSys.particles[0].get_num_springs()==1){
                num_bonded++;
            }
            else{
                num_unbonded++;
            }
            theSys.time++;
        }
        double ratio = (1.0*num_bonded)/(1.0*num_unbonded);
        std::cout << "Bonded-to-unbonded ratio:" << std::endl;
        std::cout << "Measured\tPredicted" << std::endl;
        std::cout << ratio << "\t" << std::exp(-bond_energy/theSys.kT) << std::endl;
        double prediction = std::exp(-bond_energy/theSys.kT);
        REQUIRE(std::fabs(ratio-prediction)/prediction<1e-2);
    }
}