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
        System theSys(defaultParams, myGen);
        Solver solver(theSys, defaultParams, myGen);
        REQUIRE(theSys.get_dist(theSys.particles[0],theSys.particles[1])==1.0);
        for(int t=0; t<120; t++){
            std::cout << "timestep: " << t << std::endl;
            std::cout << theSys.particles[0].get_num_springs() << std::endl;
            solver.update_bonds(theSys, theSys.dt);
            //theSys.time++;
        }
    }
}