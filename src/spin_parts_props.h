#include <math.h>
#include <sstream>
#include "glossary.h"
#include "Aboria.h"
using namespace Aboria;

class Part_set_props {

    friend class Part_set;
    
    
public:
    Part_set_props(const Glossary);

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
    double e;               // strain 
    double Fmax,Fmax2;      // Maximum force
    vdouble3 corner_0,corner_1; // Corners of box

    // @todo ; with keylist
    int init_shape;
    int init_number;
    double init_radius;
    
    bool elastic;
    double k_elast;
    
    
};
