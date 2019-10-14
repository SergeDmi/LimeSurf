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
#include "elastic_parts_props.h"
//#include "yaml-cpp/yaml.h"  // IWYU pragma: keep



Elastic_set_props::Elastic_set_props(const YAML::const_iterator config, Part_set_props & propal) : Part_set_props(propal) {
    
    init();

    Read_config(config);
}

Elastic_set_props::Elastic_set_props(const YAML::const_iterator config) : Part_set_props(config) {
    
    init();

    Read_config(config);
}

void Elastic_set_props::init() {
    prestrain=1.0;
    power_law=2;
}


void Elastic_set_props::Read_config(const YAML::const_iterator config) {
     auto conf=config->second;
    
    if (conf["prestrain"]) {
        prestrain=conf["prestrain"].as<std::double_t>();
    }
    
    if (conf["power_law"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        //power_law=conf["power_law"].as<std::int8_t>();
        power_law=conf["power_law"].as<std::double_t>();
        //power_law=2;
    }
    std::cout << "power law : " << power_law << std::endl;
    
}