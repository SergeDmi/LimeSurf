#include <math.h>
#include <sstream>
#include "glossary.h"
class Part_set_props {

    friend class Part_set;
    
    
public:
    Part_set_props(const Glossary);

    double R0;
    double minR;
    double Rmax;
    double pressure;
    double L;
    // @todo ; with keylist
    int init_shape;
    int init_number;
    
    
};