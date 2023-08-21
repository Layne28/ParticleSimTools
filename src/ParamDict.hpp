//A ParamDict consists of key-value pairs, where the keys are parameters to the simulation

#ifndef PARAM_DICT_HPP
#define PARAM_DICT_HPP

#include <string>
#include <map>

//TODO: add functionality for interpreting headers to extract subsections of conf file
//TODO: add functionality ignoring comments marked by "//"
class ParamDict
{
private:

    std::map<std::string, std::string> theMap;

public:

    //constructor
    ParamDict();

    //destructor
    ~ParamDict();

    //methods
    void read_params(std::string filename);
    void add_entry(std::string a, std::string b);
    std::string get_value(std::string key);
    bool is_key(std::string key);
    int get_size();

    //'<<' overload for std::cout
    friend std::ostream& operator<<(std::ostream& os, ParamDict& dict);
};

#endif