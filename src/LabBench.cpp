#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams, theGen), obs(theParams), solver(sys, theParams, theGen)
{

    params = theParams;
    sys.set_obs(obs);

    if(theParams.is_key("randomization_steps")) randomization_steps = std::stoi(theParams.get_value("randomization_steps"));
    if(theParams.is_key("equil_steps")) equil_steps = std::stoi(theParams.get_value("equil_steps"));
    if(theParams.is_key("production_steps")) production_steps = std::stoi(theParams.get_value("production_steps"));
    if(theParams.is_key("info_freq")) info_freq = std::stoi(theParams.get_value("info_freq"));
    if(theParams.is_key("experiment")) experiment = theParams.get_value("experiment");

    if(sys.is_aoup==1 && std::stof(theParams.get_value("va"))>0.0){
        std::cout << "WARNING: Using AOUP and active noise in the same simulation!" << std::endl;
    }
}

LabBench::~LabBench() {}

void LabBench::run(int nstps, std::string subdir, int net_freq, int therm_freq, int noisegen_freq)
{
    if (nstps==-1) nstps = this->production_steps;
    if (net_freq==-1) net_freq = this->obs.particles_freq;
    if (therm_freq==-1) therm_freq = this->obs.thermo_freq;
    if (noisegen_freq==-1) noisegen_freq = this->obs.noise_freq;

    if (obs.do_h5md==1) {
        obs.open_h5md(sys, subdir);
    }
    if (obs.do_output_noise==1) {
        obs.open_h5angen(solver, subdir);
    }

    for (int i=0; i<nstps; i++) {
        //Record data
        if (i%info_freq==0) std::cout << "step " << i << std::endl;

        if (obs.do_h5md==1) {
            if (i%net_freq==0) {
                obs.dump_h5md(sys, subdir);
            }
        }
        if (obs.do_output_noise==1) {
            if (i%noisegen_freq==0) {
                obs.dump_h5angen(solver, sys, subdir);
            }
        }
        //Advance dynamics
        solver.update(sys, sys.dt);
    }
}

void LabBench::do_experiment(std::string expt)
{
    if (expt=="standard") {
        std::cout << "Running standard experiment." << std::endl;
        this->run_standard_experiment();
    } 
    else if (expt=="ffs") {
        std::cout << "Doing forward flux sampling." << std::endl;
        this->run_ffs_experiment();
    } 
    else {
        std::cout << "This experiment has not been designed yet.\n" << std::endl;
        exit(0);
    }
}

void LabBench::run_standard_experiment()
{
    if(do_energy_minimize==1){
        std::cout << "Minimizing energy..." << std::endl;
        sys.minimize_energy();
    }

    if (randomization_steps>0) {
        std::cout << "Randomizing positions..." << std::endl;
        double kT0 = sys.kT;
        double va0 = solver.va;
        double dt0 = sys.dt;
        std::string pot = sys.nonbonded_potential_type;
        
        //Set "high" temperature and no active noise
        //to get random initial configuration
        sys.kT = 0.5;
        solver.va = 0.0;
        sys.dt = 0.00025;
        sys.nonbonded_potential_type = "wca";
        std::cout << "kT: " << sys.kT << std::endl;;
        std::cout << "va: " << solver.va << std::endl;
        std::cout << "potential: " << sys.nonbonded_potential_type << std::endl;
        std::cout << "dt: " << sys.dt << std::endl;;
        this->run(this->randomization_steps, "/randomize", this->obs.particles_freq, this->obs.thermo_freq);

        //reset parameters for production
        sys.kT = kT0;
        solver.va = va0;
        sys.dt = dt0;
        sys.nonbonded_potential_type = pot;

        std::cout << "Resetting parameters: " << std::endl;
        std::cout << "kT: " << sys.kT << std::endl;
        std::cout << "va: " << solver.va << std::endl;
        std::cout << "potential: " << sys.nonbonded_potential_type << std::endl;
    }
    
    std::cout << "Equilibrating..." << std::endl;
    this->run(this->equil_steps, "/equil", this->obs.particles_freq, this->obs.thermo_freq);

    std::cout << "Doing production run..." << std::endl;
    this->run(this->production_steps, "/prod", this->obs.particles_freq, this->obs.thermo_freq);
}

auto LabBench::run_ffs_stage1(int N0, double la, std::string op)
{
    //This function has to be placed before "run_ffs_experiment"
    //because of the use of auto return type

    //Define struct for output of stage 1
    struct result {
        double time; //total time to collect N0 configurations
        std::vector<std::vector<Particle>> configs; //vector of configurations crossing la
    };

    std::vector<std::vector<Particle>> configs;
    configs.push_back(sys.particles);

    return result {1.0, configs};
}

void LabBench::run_ffs_experiment()
{
    //Do forward flux sampling (ffs) using original "direct" algorithm
    //INPUT:
    //  -N0 (number of points at first interface)
    //  -M0 (number of "firing runs" at first interface)
    //  -nint (number of interfaces)
    //  -la (op value demarking edge of state A)
    //  -lb (op value demarking edge of state B)
    //  -op (order parameter)

    int N0 = 5;
    int M0 = 5;
    int nint = 5;
    double la = 0.0;
    double lb = 1.0;
    std::string op = "single_particle_x";

    if(params.is_key("N0")) N0 = std::stoi(params.get_value("N0"));
    if(params.is_key("M0")) M0 = std::stoi(params.get_value("M0"));
    if(params.is_key("nint")) nint = std::stoi(params.get_value("nint"));
    if(params.is_key("la")) la = std::stod(params.get_value("la"));
    if(params.is_key("lb")) lb = std::stod(params.get_value("lb"));
    if(params.is_key("op")) op = params.get_value("op");

    //Stage 1: sampling the A basin

    //TODO: reset system to make sure it starts in A basin
    auto [Tstage1, stage1_configs] = run_ffs_stage1(N0, la, op);

    std::cout << "Tstage1: " << Tstage1 << std::endl;

    //Stage 2: crossing interfaces

}

    /*
    std::cout << "Testing grid" << std::endl;
    
    std::array<double,2> min = {-0.5*sys.edges[0], -0.5*sys.edges[1]};
    std::array<double,2> max = {0.5*sys.edges[0], 0.5*sys.edges[1]};
    std::array<bool,2> per = {true, true};
    NeighborGrid<Particle, 2> *grid = new NeighborGrid<Particle,2>(min, max, per, 3.5);
    //grid->update_atom(&sys.particles[0]);
    //grid->update_atom(&sys.particles[1]);
    //grid->update_atom(&sys.particles[2]);
    //auto p = grid->get_neighbor_iterator(&sys.particles[0]);
    //std::cout << *p << std::endl;
    Particle p1;
    p1.pos[0]=0.0;
    p1.pos[1]=0.0;
    Particle p2;
    p2.pos[0] = 1.0;
    p2.pos[1] = 0.0;
    grid->update_atom(&p1);
    grid->update_atom(&p2);
    for (auto p = grid->get_neighbor_iterator(&p1); !p.isDone(); p++){ 
        std::cout << *p << std::endl;
    }
    std::cout << "passed" << std::endl;
    */