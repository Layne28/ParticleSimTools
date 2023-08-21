#include <cstdlib>
#include <iostream>
#include <fstream>
#include <experimental/filesystem>
#include <cmath>

#include "LabBench.hpp"

using namespace std;
namespace fs = experimental::filesystem;
void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams);

int main(int argc, char * argv[])
{
    //defaults
    std::string input_file = "sample.in";
    std::string seed_file = "";
    long unsigned int seed = 1;
    if (argc>1) input_file = argv[1];
    if (argc>2) seed = std::atoi(argv[2]);
    if (argc>3) seed_file = argv[3];
    std::cout << "Using input file: " << input_file << std::endl;

    ParamDict myParams;
    myParams.read_params(input_file);

    //Seed RNG
    determine_seed(seed, seed_file, myParams);
    std::cout << myParams.get_value("output_dir") << std::endl;

    std::cout << "Using seed: " << seed << std::endl;

    gsl_rng *myGen = CustomRandom::init_rng(seed);

    LabBench myBench(myParams, myGen);
    myBench.do_experiment(myBench.experiment);

    return 0;
}

void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams)
{
    if(seed_file!="" && fs::exists(seed_file)){
        myParams.add_entry("output_dir", myParams.get_value("output_dir") + "/seed=" + std::to_string(seed));
        std::cout << "Reading seed from file." << std::endl;
        std::ifstream file(seed_file);
        file.exceptions(std::ifstream::failbit);
        int found = 0;
        try
        {
            std::string line;
            int cnt = 1;
            while(std::getline(file,line))
            {
                int new_seed = std::stoi(line);
                if(seed==(long unsigned)cnt){
                    seed = (long unsigned) new_seed;
                    found = 1;
                    break;
                }
                cnt += 1;
            }

        }
        catch (std::ifstream::failure &e)
        {
            if (!file.eof()) std::cout << "Exception opening/reading/closing seed file.\n";
        }
        if(found==0){
            std::cout << "Error: could not find seed in file!\n";
            exit(-1);
        }
    }
    fs::create_directories(myParams.get_value("output_dir"));
    std::ofstream ofile;
    ofile.open(myParams.get_value("output_dir") + "/seed_value.txt");
    ofile << seed << std::endl;
    ofile.close();
}