#include <math.h>
#include <sstream>
#include "glossary.h"

#ifndef MESHLESS_SPIN_PROPS_H
#define MESHLESS_SPIN_PROPS_H


class Meshless_props {
    
public:
    // constructor
    Meshless_props(const Glossary);
    
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
    };
    
    // Time sequence
    typedef std::vector<double> double_list;
    double_list change_times;
    double_list change_press;
    int nb_times;
    
};

#endif /* MESHLESS_SPIN_PROPS_H */

