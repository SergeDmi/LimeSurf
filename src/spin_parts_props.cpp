#include <random>
//#include "glossary.h"
//#include "Aboria.h"
//using namespace Aboria;
#include <sstream>
#include "spin_parts_props.h"
using namespace std;

// Empty constructor for convenience
Part_set_props::Part_set_props() {
    init();
};


/*
Part_set_props::Part_set_props(const Glossary opt) {
    // @TODO : replace all glos - YAML
    init();
    
    opt.set(renorm_rate, "renorm_rate");
 
    
    // Confinement 
    //x_conf=0;y_conf=0;z_conf=0;
    //x_max=0 ;y_max=0 ;z_max=0 ;
    
    // reading these parameters from the config file !
    /*
    std::string str;
    str.reserve(10000);
    fname_in.reserve(1000);
    //load_from_file=opt.set(str,"load_from");
    fname_in.append(str);
        //std::cout << "fname_in : " << fname_in << std::endl;
    
    if (load_from_file>0) {
        fname_in.append(str);
        std::cout << "fname_in : " << fname_in << std::endl;
    }
    
    opt.set(z_max, "z_conf",0);
    opt.set(z_conf, "z_conf",1);
    
    opt.set(y_max, "y_conf",0);
    opt.set(y_conf, "y_conf",1);
    
    opt.set(x_max, "x_conf",0);
    opt.set(x_conf, "x_conf",1);
   
    
    opt.set(mechanics, "mechanics");
    
    opt.set(k_elast, "k_elast");
    
}

*/

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
        std::cout << " set visco to " << visco << std::endl;
    }
    
     if (conf["Rvisc"]) {
        Rvisc=conf["Rvisc"].as<std::double_t>();
    }
    
    if (conf["type"]) {
        //pressure=conf["pressure"].as<std::double_t>();
        mechanics=conf["type"].as<std::double_t>();
    }
 
}