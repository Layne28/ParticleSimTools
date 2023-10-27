#include "catch.hpp"
#include "../../src/Solver.hpp"
#include "../../src/Particle.hpp"
#include "../../src/Spring.hpp"
#include "../../src/System.hpp"
#include "../../src/CustomRandom.hpp"

TEST_CASE("Solver")
{
    ParamDict defaultParams;
    gsl_rng *myGen = CustomRandom::init_rng(1);
    defaultParams.add_entry("k0_bond", "1.0");
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
        double dt;
        double lifetimes[nsim] = {0};
        for(int s=0; s<nsim; s++){
            System theSys(defaultParams, myGen);
            Solver solver(theSys, defaultParams, myGen);
            dt = theSys.dt;

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
        //for(int s=0; s<nsim; s++) std::cout << "Lifetime: " << lifetimes[s] << std::endl;

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
        //double testsum = 0.0;
        //double testsum2 = 0.0;
        for(int i=0; i<nbins/3; i++){
            std::cout << i*binwidth << "\t" << histogram[i] << "\t" << std::exp(-(i+0.5)*binwidth)*binwidth << std::endl;
            //Check for more than than 10% relative error in bins with significant probability
            //If this fails, try increasing "nsim"
            if(histogram[i]>0.3*histogram[0]){
                REQUIRE(std::fabs(histogram[i]-std::exp(-(i+0.5)*binwidth)*binwidth)/(std::exp(-(i+0.5)*binwidth)*binwidth)<0.20); 
            }
            //testsum += std::exp(-(i)*binwidth)*binwidth;
            //testsum2 += histogram[i];
        }
        //std::cout << "test: " << testsum << std::endl;
        //std::cout << "test2: " << testsum2 << std::endl;
    }
}