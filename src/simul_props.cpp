#include <sstream>
#include "simul_props.h"
#include "Aboria.h"
using namespace Aboria;

Simul_props::Simul_props() {
    init();
}

Simul_props::Simul_props(const YAML::const_iterator config) {
  
    init();
    vdouble3 zero(0,0,0);

  
    Read_config(config);

    // Number of times we should save the shape
    
    
    pressure=0;
    confine_pot=zero;
    confine_max=zero;
    confine_min=zero;
    //x_conf=0;y_conf=0;z_conf=0;
    
    Read_config(config);
}

void Simul_props::init() {
    // Times
    dt=0.01;
    Tend=1.0;
    n_frames=2;
    name="";
    nb_times=0;
    
    // By default, no confinement of course
    is_confinement=0;
    add_confinement_force=[] (vdouble3& force, vdouble3 posi ) {};
}

void Simul_props::Read_config(const YAML::const_iterator config) {
    //auto conf=config->second;
    //opt.set(dt,"dt");
    //opt.set(Tend,"Tend");
    //opt.set(n_frames,"n_frames");
    name=config->first.as<std::string>();
    auto conf=config->second;
    
    if (conf["duration"]) {
        std::cout << "getting duration from : " << conf["duration"].as<std::double_t>() << std::endl;
        Tend=conf["duration"].as<std::double_t>();
    }
    
     if (conf["dt"]) {
         std::cout << "getting dt from : " << conf["dt"].as<std::double_t>() << std::endl;
        dt=conf["dt"].as<std::double_t>();
    }
    
     if (conf["n_frames"]) {
         std::cout << "getting n_frames from : " << conf["n_frames"].as<std::double_t>() << std::endl;
        n_frames=conf["n_frames"].as<std::double_t>();
    }
    
    if (conf["pressure"]) {
        pressure=conf["pressure"].as<std::double_t>();
    }
    if (conf["confinement"]) {
        int i=-1;
        YAML::Node confs=conf["confinement"];
        for(auto confine=confs.begin();confine!=confs.end();++confine) {
            std::string key=confine->first.as<std::string>();
            // std::cout << "WESH we got " << key << std::endl;
            if        (key=="x") {
                i=0;
            } else if (key=="y") {
                i=1;
            } else if (key=="z") {
                i=2;
            }
            
             if (i>=0) {
                //std::cout << "somehow got that we are confied on axis : " << i << std::endl;
                YAML::Node coco=confine->second;
                if        (coco["min"]) {
                    confine_min[i]=coco["min"].as<std::double_t>();
                }
                if (coco["max"]) {
                    confine_max[i]=coco["max"].as<std::double_t>();
                }
                if (coco["stiffness"]) {
                    confine_pot[i]=coco["stiffness"].as<std::double_t>();
                    //std::cout << key << " : " << confine_min[i] << " , " << confine_max[i] << " : " << confine_pot[i] << std::endl;
                }
            }   
            else {
                std::cout << "Warning : could not understand confinement axis" << key << std::endl;
            }

        }
       
            
        //    confine_pot[i]=confs[i].as<std::double_t>();
        //pressure=conf["confinements"].as<std::double_t>();
    }
    
    
    
    if (n_frames) {
        dt_frames=Tend/(n_frames);
    }
    else {
        dt_frames=Tend;
    }
    n_frames++;
   
    Make_confinement();
}




void Simul_props::Make_confinement() {
    std::vector <int>is_confined;
    int any_conf=0;
    for (int i=0;i<3;++i) {
        if (abs(confine_pot[i])>0) {
            any_conf++;
            is_confined.push_back(i);
        }
    }
    
    if (any_conf>0) {
        is_confinement=1;
        
        add_confinement_force=
         [conf=confine_pot,xmax=confine_max,xmin=confine_min,ixes=is_confined] (vdouble3& force, vdouble3 posi ) {
            for (auto j=ixes.begin();j!=ixes.end();++j) {
                force[*j]-=conf[*j]*( 
                              ( posi[*j] > xmax[*j] )*( posi[*j] - xmax[*j] ) 
                            + ( posi[*j] < xmin[*j] )*( posi[*j] - xmin[*j] )  
                        ); 
            }
         };
    }
    
}