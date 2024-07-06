#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams, theGen), obs(theParams), solver(sys, theParams, theGen)
{

    params = theParams;
    rg = theGen;
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

/******************************/
/*** Forward Flux Sampling ****/
/******************************/

auto LabBench::run_ffs_stage1(int N0, double l0, double la, double lb)
{
    //This function has to be placed before "run_ffs_experiment"
    //because of the use of auto return type
    //TODO: handle case where system reaches state B (rare but possible)

    //Define struct for output of stage 1
    struct result {
        double time; //total time to collect N0 configurations
        std::vector<std::vector<Particle>> configs; //vector of configurations crossing la
    };

    //Define variables
    std::vector<std::vector<Particle>> configs; //vector of configurations

    int config_counter = 0; //no. of configs reaching la (out of N0)
    double time = 0;        //total time to get N0 configs

    int was_in_a = 1;       //keep track of whether system was in state A
                            //before crossing l0

    //Dynamics loop
    int nsteps = 0;
    int outfreq = 100000;
    while (config_counter<N0){

        if (nsteps % outfreq == 0){
            std::cout << "time: " << time << " timesteps: " << nsteps << std::endl;
        }

        //Run one step
        solver.update(sys, sys.dt);
        time += sys.dt;
        nsteps++;

        if (sys.get_order_parameter() <= la){
            was_in_a = 1;
        }
        //Append configuration and increment counter if 
        //system has crossed lambda_a
        if (sys.get_order_parameter() > l0 && was_in_a==1){
            std::cout << "particle crossed lambda_0 at time " << time << std::endl;
            configs.push_back(sys.particles);
            config_counter++;
            was_in_a = 0;
        }
    }
    

    return result {time, configs};
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
    if(params.is_key("order_parameter")) op = params.get_value("order_parameter");

    sys.order_parameter = op;

    //Define interfaces
    //For now, use simplest method: nint evenly spaced
    //interfaces between la and lb
    std::vector<double> lambdas(nint, 0.0);
    for(int i=0; i<nint; i++){
        lambdas[i] = la + (i+1)*(lb-la)/nint;
        std::cout << lambdas[i] << std::endl;
    }

    /*****************************/
    //Stage 1: sampling the A basin
    /*****************************/

    //TODO: reset system to make sure it starts in A basin
    auto [Tstage1, stage1_configs] = run_ffs_stage1(N0, lambdas[0], la, lb);

    //Compute flux from A to 0
    double phi_A0 = N0/Tstage1;
    std::cout << "Phi_A,0: " << phi_A0 << std::endl;

    /*****************************/
    //Stage 2: crossing interfaces
    /*****************************/
    
    //TODO: save configurations of successful trajectories 

    std::vector<std::vector<Particle>> stage2_configs_prev = stage1_configs;
    std::vector<double> transition_probs(nint-1, 0);

    //Loop through interfaces
    for(int i=1; i<nint; i++){

        std::cout << "transition from lambda_" << (i-1) << " to lambda_" << i << std::endl;

        std::vector<std::vector<Particle>> stage2_configs_curr;
        int Ni=0; //number of trajectories crossing lambda_(i+1)

        //Launch M_i (=M0 for now) trial trajectories
        //randomly sampled from the N_i (=N0 for now) stored configurations
        for(int j=0; j<M0; j++){

            std::cout << "test" << std::endl;
            int index = gsl_rng_uniform_int(rg, N0);
            std::cout << index << std::endl;
            std::vector<Particle> config_curr = stage2_configs_prev[index];

            //Assign system the selected configuration
            //TODO: if pbc, may also need to store and select
            //periodic images, update cell list, etc
            sys.particles = config_curr;

            //Dynamics loop
            int done = 0;
            while (done==0){
                solver.update(sys, sys.dt);

                if (sys.get_order_parameter() <= la){
                    std::cout << "system returned to state A" << std::endl;
                    done = 1;
                }
                
                if (sys.get_order_parameter() > lambdas[i]){
                    std::cout << "particle crossed lambda_" << i << std::endl;
                    stage2_configs_curr.push_back(sys.particles);
                    done = 1;
                    Ni++;
                }
            }
        }
        transition_probs[i-1] = (1.0*Ni)/M0;
        stage2_configs_prev = stage2_configs_curr;
    }

    std::cout << "Transition probabilities:" << std::endl;
    for(int i=0; i<(nint-1); i++){
        std::cout << i << "-->" << (i+1) << ": " << transition_probs[i] << std::endl;
    }

    double kAB = phi_A0;
    for (int i=0; i<(nint-1); i++){
        kAB *= transition_probs[i];
    }
    std::cout << "Rate constant: " << kAB << std::endl;



}

/******************************/
/*** Umbrella Sampling ********/
/******************************/