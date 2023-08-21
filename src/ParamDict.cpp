#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include "ParamDict.hpp"

ParamDict::ParamDict() {}

ParamDict::~ParamDict() {}

int ParamDict::get_size()
{
    return this->theMap.size();
}

std::string ParamDict::get_value(std::string key)
{
    return this->theMap[key];
}

bool ParamDict::is_key(std::string key)
{
    return this->theMap.count(key);
}

void ParamDict::add_entry(std::string a, std::string b)
{
    this->theMap[a] = b;
}

void ParamDict::read_params(std::string filename)
{
    //Check file extension
    int delim_index = filename.find_last_of(".");
    std::string extension = filename.substr(delim_index+1);
    if (extension!="in")
    {
        throw std::runtime_error("Error: parameters must be read from a .in file!");
    }
    std::ifstream file(filename);
    file.exceptions(std::ifstream::failbit);// | std::ifstream::badbit);
    try
    {
        //Read key-value pairs in line-by-line
        std::string line;
        while(std::getline(file,line))
        {
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end()); //remove whitespace
            int index = line.find("=");
            if (index==-1) continue;

            std::string key = line.substr(0,index);
            std::string value = line.substr(index+1);
            this->add_entry(key, value);
            
        }
    }
    catch (std::ifstream::failure &e)
    {
        if (!file.eof()) std::cerr << "Exception opening/reading/closing in file\n";
    }
    
}

std::ostream& operator<<(std::ostream& os, ParamDict& dict)
{
    for (const auto &[k,v] : dict.theMap)
    {
        os << k << "=" << v << std::endl;
    }
    return os;
}