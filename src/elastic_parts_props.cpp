/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include <sstream>
#include "elastic_parts_props.h"




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
    young_modulus=-1;
    thickness=0.0;
}


void Elastic_set_props::Read_config(const YAML::const_iterator config) {
     auto conf=config->second;
    
    if (conf["prestrain"]) {
        prestrain=conf["prestrain"].as<std::double_t>();
    }
    
    if (conf["power_law"]) {
        power_law=conf["power_law"].as<std::double_t>();
    }
    
    if (conf["young_modulus"]) {
        young_modulus=conf["young_modulus"].as<std::double_t>();
    }
    
    if (conf["thickness"]) {
        thickness=conf["thickness"].as<std::double_t>();
    }
        
}


bool Elastic_set_props::Check_props() {
    return Part_set_props::Check_props();
}
