#include <math.h>
#include <sstream>
#include "glossary.h"
class Part_set_props {

    friend class Part_set;
    
    
public:
    Part_set_props(const Glossary);

    double R0;
    double minR;
    double minR2;
    double Rmax;
    double pressure;
    double k_rep;
    double k_att;
    double p_rep;
    double p_att;
    double k_bend;
    double k_align;
    double visco;
    double Rvisc;
    double L;
    double Fmax;
    // @todo ; with keylist
    int init_shape;
    int init_number;
    double init_radius;
    
    
};
