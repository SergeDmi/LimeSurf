/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include <sstream>
#include "spin_parts_props.h"
using namespace std;

// Empty constructor for convenience
Part_set_props::Part_set_props() {
    init();
};

Part_set_props::Part_set_props(const YAML::const_iterator config) {
    
    init();
    Read_config(config);
}

Part_set_props::Part_set_props(const Part_set_props & props) {
    
    init();
    fname_out=props.fname_out;
    fname_in=props.fname_in;
    mechanics=props.mechanics;
    k_align=props.k_align;
    k_elast=props.k_elast;
    visco=props.visco;
    Rvisc=props.Rvisc;
    renorm_rate=props.renorm_rate;
    verbose=props.verbose;
}

void Part_set_props::init() {
      
    fname_out="simulated_";
    fname_in="config.yaml";
    // Shape to be created
    
    mechanics=1;
    verbose=0;
    
    // Physical parameters
    // For any set :
    k_align=1.0;
    k_elast=1.0;
    visco=1.0;
    Rvisc=1.0;
	mobility=0.0;
    //pressure=0;
    //relax=0;
    
   
    
    // Misc parameters
    renorm_rate=0.001;
}


//void Part_set_props::Read_config(const YAML::const_iterator config) {
    
    
void Part_set_props::Read_config(const YAML::const_iterator config) {
    name=config->first.as<std::string>();
    fname_out=name;
    
    auto conf=config->second;
    
    if (conf["source"]) {
        std::string str;
        str.reserve(1000);
        fname_in.reserve(1000);
        str=conf["source"].as<std::string>();
        fname_in.assign(str);
        
    }
    
    if (conf["out"]) {
        std::string str1;
        str1.reserve(100);
        fname_in.reserve(100);
        str1=conf["out"].as<std::string>();
        fname_out.assign(str1);
    }
    
    if (conf["k_elast"]) {
        k_elast=conf["k_elast"].as<std::double_t>();
    }
    
    
    if (conf["align"]) {
        k_align=conf["align"].as<std::double_t>();
    }
    
     if (conf["renorm"]) {
        renorm_rate=conf["renorm"].as<std::double_t>();
    }
    
     if (conf["verbose"]) {
        verbose=conf["verbose"].as<std::double_t>();
    }
    
	if (conf["visco"]) {
        visco=conf["visco"].as<std::double_t>();
        //std::cout << " set visco to " << visco << std::endl;
    }
	
	if (conf["viscosity"]) {
        visco=conf["viscosity"].as<std::double_t>();
        //std::cout << " set visco to " << visco << std::endl;
    }
	
	if (conf["mobility"]) {
        mobility=conf["mobility"].as<std::double_t>();
        //std::cout << " set visco to " << visco << std::endl;
    }
    
    if (visco<0) {
        std::cerr << "Setting viscosity to infinity" << std::endl;
        visco=INFINITY;
    }
    
     if (conf["Rvisc"]) {
        Rvisc=conf["Rvisc"].as<std::double_t>();
    }
    
    
    if (Rvisc<0) {
        std::cerr << "Setting rot. viscosity to infinity" << std::endl;
        Rvisc=INFINITY;
    }
    
    if (conf["type"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        mechanics=conf["type"].as<std::double_t>();
    }
 
}

void Part_set_props::PrintIfVerbose(const std::string & message) const {
    if (verbose) {
        std::cout << "[ " << name << " ] : " << message << std::endl;
    }
}