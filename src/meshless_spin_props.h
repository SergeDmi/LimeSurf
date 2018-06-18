#include <math.h>
#include <sstream>
#include "glossary.h"

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
};
