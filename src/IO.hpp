#ifndef IO_HPP
#define IO_HPP

#include "ParamDict.hpp"
#include <iostream>
#include <fstream>
#include <experimental/filesystem>

void determine_seed(long unsigned int &seed, std::string seed_file, ParamDict &myParams);
#endif