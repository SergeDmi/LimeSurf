#include <math.h>
#include <sstream>
#include "glossary.h"
#include "Aboria.h"
using namespace Aboria;
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

#ifndef SIMUL_PROPS_H
#define SIMUL_PROPS_H


class Simul_props {
    
public:
    // Dummy constructor
    Simul_props();

    // constructor
    // Deprecated : Simul_props(const Glossary) ;
    Simul_props(const YAML::const_iterator);
    
    
    // Initialization
    void init();
    
    // Read config
    //void Read_config(const YAML::Node,const YAML::Node);
    void Read_config(const YAML::const_iterator);

    
    // Update 
    void Make_confinement();
    
    // Name
    std::string name;
    
    // Output name
    std::string out_name;
    
    // Simul time step
    double dt;
    
    // End of simulation time
    double Tend;
    
    // Duration of a frame (between two saves)
    double  dt_frames;
    
    // Number of frames (saving timepoints)
    int n_frames;
  
    // Type of simulation
    int mechanics;    
    
    enum SimulType
    {
        VISCOUS  = 0,           // Purely viscous simulation
        ELASTIC  = 1,           // Purely elastic simulation
        SIMPLE_VISCOEL = 2,     // Simple visco-elastic implementation
        TETRA_ELASTIC  = 4,     // Purely elastic tetrahedron simulation
    };
    
    // Time sequence
    typedef std::vector<double> double_list;
    double_list change_times;
    double_list change_press;
    int nb_times;
    
    // Pressure
    double pressure; 
    
    // Confinement
    
     // Confinement forces
    std::function<void(vdouble3& ,const vdouble3)> add_x_conf;
    std::function<void(vdouble3& ,const vdouble3)> add_y_conf;
    std::function<void(vdouble3& ,const vdouble3)> add_z_conf;
    std::function<void(vdouble3& ,const vdouble3)> add_confinement_force;
    //double conf_xmax,conf_xmin,conf_ymax,conf_ymin , conf_zmax,conf_zmin;
    //double x_conf,y_conf,z_conf;
    bool is_confinement;
    
    
private:
    vdouble3 confine_pot;
    vdouble3 confine_max;
    vdouble3 confine_min;
    
    
};

#endif /* MESHLESS_SPIN_PROPS_H */

