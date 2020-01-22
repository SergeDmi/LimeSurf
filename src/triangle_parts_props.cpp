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
#include "triangle_parts_props.h"

/*
Triangle_set_props::Triangle_set_props(const Glossary opt) : Elastic_set_props(opt) {
    
    prestrain=1.0;
    
    opt.set(prestrain, "prestrain");
    
    rigidity=0.0;
    
    opt.set(rigidity, "flexural");
}
*/



Triangle_set_props::Triangle_set_props(const YAML::const_iterator config) : Elastic_set_props(config) {
    
    auto conf=config->second;
    
    k_bending=-1;
    
    imposed_angle=-1;
    
    if (conf["bending"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        k_bending=conf["bending"].as<std::double_t>();
    }
    
    if (conf["k_bending"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        k_bending=conf["k_bending"].as<std::double_t>();
    }
    
    if (conf["imposed_angle"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        imposed_angle=conf["imposed_angle"].as<std::double_t>();
    }
    std::cout << "bending : " << k_bending << std::endl;
    
}

