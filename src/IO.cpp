#include "IO.hpp"

using namespace std;
namespace fs = experimental::filesystem;

void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams)
{
    myParams.add_entry("output_dir", myParams.get_value("output_dir") + "/seed=" + std::to_string(seed));
    if(seed_file!="" && fs::exists(seed_file)){
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