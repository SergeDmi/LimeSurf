#include <math.h>
#include <sstream>
#include "glossary.h"
#include "Aboria.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace Aboria;

#ifndef SPIN_PARTS_PROPS_H
#define SPIN_PARTS_PROPS_H

class Part_set_props {

    friend class Part_set;
    friend class Elastic_part_set;
    friend class Simple_viscoel_part_set;
    
public:
    // Empty constructor
    Part_set_props();
    
    // Real constructor
    Part_set_props(const Glossary);
    std::string fname_in;   // contains file to be loaded
    std::string fname_out;
    //fname_in.reserve(10000);
    int load_from_file;    // flag if a file needs to be loaded
    double R0;              // Equilibrium distance between parts
    double minR,minR2;      // Minimum interaction distance (avoid dvg)
    double Rmax;            // Max interaction distance
    double Rsearch;         // Search radius
    double pressure;        // A pressure
    double k_rep,p_rep;     // repulsive potential and power law
    double k_att,p_att;     // attractive potential and power law
    double k_bend,p_bend;   // bending potential and power law
    double k_align,p_align; // alignment potential and power law
    double visco,Rvisc;     // viscosity (translational & rotational)
    double L;               // box size
    
    double renorm_rate;     // rate at which normals are normalized can be quite low
    double Fmax,Fmax2;      // Maximum force
    double relax;           // linkers relaxation rate
    vdouble3 corner_0,corner_1; // Corners of box
    int max_neighbours;
    
    double x_max, y_max, z_max ;
    double x_conf,y_conf,z_conf;
    
    // @todo ; with keylist
    int init_shape;
    int init_number;
    double init_radius;
    
    bool elastic;
    double k_elast;
    
    virtual void read_physical_properties(YAML::Node);
    
};

#endif