/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/

#include <random>
#include <sstream>
#include "simple_viscoel_parts_props.h"

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