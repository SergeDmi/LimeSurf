#include <math.h>
#include <sstream>
//#include "glossary.h"
#include "Aboria.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace Aboria;

#ifndef SPIN_PARTS_PROPS_H
#define SPIN_PARTS_PROPS_H

// @TODO :
// implement copy constructor ?

class Part_set_props {

    friend class Part_set;
    friend class Elastic_part_set;
    friend class Simple_viscoel_part_set;
    
public:
    // Dummy constructor
    Part_set_props();
    
     // constructor
    Part_set_props(const YAML::const_iterator);
    
    // Copy constructor
    Part_set_props(const Part_set_props &);
    
    // Initialization
    void init();
    
    // Actual reading from file
    void Read_config(const YAML::const_iterator);
    
    // Type of simulation
    int mechanics;    
    
    enum MechType
    {
        VISCOUS  = 0,           // Purely viscous simulation <- NOT IMPLEMENTED
        ELASTIC  = 1,           // Purely elastic simulation
        SIMPLE_VISCOELASTIC = 2,     // Simple visco-elastic implementation <- NOT IMPLEMENTED
        TETRA_ELASTIC  = 4,     // Purely elastic tetrahedron simulation
    };

    std::string fname_in;   // contains file to be loaded
    std::string fname_out;
    std::string name;
    
    //int load_from_file;    // flag if a file needs to be loaded
    double k_elast ;            // Spring constant
    double k_align ;            // normal alignment constant
    double renorm_rate;         // Rate at which we renormalize vectors (higher is better but slower)
    double visco ;              // translational viscosity
    double Rvisc ;              // rotational viscosity of normals
    
    // Do we talk about it ?
    bool verbose;
    
};

#endif