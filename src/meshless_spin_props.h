#include <math.h>
#include <sstream>
#include "glossary.h"

class Meshless_props {
    
public:
    Meshless_props(const Glossary);

    double dt;
    double Tend;
    double  dt_frames;
    int n_frames;
};