/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <random>
//#include "glossary.h"
//#include "Aboria.h"
//using namespace Aboria;
#include <sstream>
#include "simple_viscoel_parts_props.h"

/*
Simple_viscoel_set_props::Simple_viscoel_set_props(const Glossary opt) : Elastic_set_props(opt) {
    
    time_plastic=DBL_MAX;
    size_plastic=DBL_MAX;
    
    opt.set(time_plastic, "plasticity",0);
    opt.set(size_plastic, "plasticity",1);
    
    
}
*/


Simple_viscoel_set_props::Simple_viscoel_set_props(const YAML::const_iterator config) : Elastic_set_props(config) {
    auto conf=config->second;
    
    time_plastic=DBL_MAX;
    size_plastic=DBL_MAX;
    
    if (conf["time_plastic"]) {
        time_plastic=conf["time_plastic"].as<std::double_t>();
    } else {
        std::cout << "WARNING : no time_plastic found" << std::endl;
    }
    
    if (conf["size_plastic"]) {
        size_plastic=conf["size_plastic"].as<std::double_t>();
    } else {
        std::cout << "WARNING : no size_plastic found" << std::endl;
    }
    
    
}