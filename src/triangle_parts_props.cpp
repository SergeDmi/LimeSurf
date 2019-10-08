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
    
   prestrain=1.0;
   rigidity=1.0;
   
    auto conf=config->second;
    
    if (conf["prestrain"]) {
        prestrain=conf["prestrain"].as<std::double_t>();
    }
    
    if (conf["flexural"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        rigidity=conf["flexural"].as<std::int8_t>();
    }
    //std::cout << "power law : " << power_law << std::endl;
    
}

