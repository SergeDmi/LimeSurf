#include <random>
//#include "glossary.h"
//#include "Aboria.h"
//using namespace Aboria;
#include <sstream>
#include "viscous_parts_props.h"
using namespace std;

// Empty constructor for convenience
Viscous_set_props::Viscous_set_props() : Part_set_props() {
    std::cerr << "ERROR : Viscous set not implemented" << std::endl;
    init();
};


/*
Viscous_set_props::Viscous_set_props(const Glossary opt) : Part_set_props(opt) {
    // @TODO : replace all glos - YAML
    init();
    
    opt.set(renorm_rate, "renorm_rate");
    double box[6];
    
    // Confinement 
    //x_conf=0;y_conf=0;z_conf=0;
    //x_max=0 ;y_max=0 ;z_max=0 ;
    
    // reading these parameters from the config file !
    
    std::string str;
    str.reserve(10000);
    fname_in.reserve(1000);
    load_from_file=opt.set(str,"load_from");
    
    if (load_from_file>0) {
        fname_in.append(str);
        std::cout << "fname_in : " << fname_in << std::endl;
    }
     
    
    opt.set(init_shape, "shape", KeyList<int>("sphere", 0, "sheet", 1, "pombe", 2));
    opt.set(init_number, "number");
    opt.set(init_radius, "radius");
    opt.set(relax, "relax");
    
    opt.set(visco, "visco",0);
    opt.set(Rvisc, "visco",1);
    opt.set(Rvisc, "Rvisco");
    opt.set(max_neighbours,"max_neighbours");
    
   
    opt.set(k_att, "attraction",0);
    opt.set(k_rep, "repulsion",0);
    opt.set(p_att, "attraction",1);
    opt.set(p_rep, "repulsion",1);
    
    p_bend=p_att;
    p_align=p_att;
    
    opt.set(k_align, "align",0);
    
    /*
    opt.set(z_max, "z_conf",0);
    opt.set(z_conf, "z_conf",1);
    
    opt.set(y_max, "y_conf",0);
    opt.set(y_conf, "y_conf",1);
    
    opt.set(x_max, "x_conf",0);
    opt.set(x_conf, "x_conf",1);
    
     
     
    opt.set(k_bend, "bending",0);
    opt.set(p_align, "align",1);
    opt.set(p_bend, "bending",1);
    
    //opt.set(elastic, "elastic");
    
    for (int i=0; i<6; ++i) {
        opt.set(box[i], "box",i);
    }
    

    corner_0=vdouble3(box[0],box[1],box[2]);
    corner_1=vdouble3(box[3],box[4],box[5]);
    
    // Computing R0
    if (k_att>0) {
            R0=pow(p_rep*k_rep/(p_att*k_att),1.0/(p_rep-p_att));
    }
    opt.set(R0, "R0");
    minR=0.75*R0;
    Rmax=2.0*R0;

    k_elast=1.0;
    opt.set(k_elast, "k_elast");
    opt.set(minR, "Rmin");
    minR2=pow(minR,2.0);
    if (k_rep>0) {
        Fmax=k_rep/pow(minR,p_rep);
    }
    opt.set(Fmax, "Fmax");
    Fmax2=Fmax*Fmax;
    
    opt.set(Rmax, "Rmax");
    Rsearch=Rmax;
    opt.set(Rsearch, "Rsearch");
    opt.set(Rmax, "Rmax");
    
    //opt.set(pressure, "pressure");
    
    std::cout << "# Fmax2=" << Fmax2 << std::endl;
    
}

 */
 
 
Viscous_set_props::Viscous_set_props(const YAML::const_iterator config) : Part_set_props(config) {
    init();
    
    Read_config(config);
}


void Viscous_set_props::init() {
    Part_set_props::init();
    fname_out="simulated_";
    // Shape to be created
    init_shape=0;
    init_number=100;
    init_radius=50;
    
    // kind of set
    //elastic=0;
    
    // Distance between particles & box size
    // @TODO : this needs a serous cleanup
    Rmax=3;
    minR=0.5;
    L=200.0;
    
    // Physical parameters
    // For any set :
    k_align=1.0;
    k_bend=1.0;
    visco=0.0;
    viscosity=1.0
   
    //pressure=0;
    relax=0;
    
    // For viscous sets
    k_att=2.0;
    k_rep=1.0;
    p_rep=12.0;
    p_att=6.0;
    R0=1.0;
    Fmax=DBL_MAX;
    
    // For elastic sets :
    max_neighbours=15;
    
    // Misc parameters
    renorm_rate=0.001;
}


//void Part_set_props::Read_config(const YAML::const_iterator config) {
    
    
void Viscous_set_props::Read_config(const YAML::const_iterator config) {
    Part_set_props::Read_config(config);
}