#include <random>
#include "Aboria.h"
using namespace Aboria;
#include "spin_parts.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>


#include "spin_parts_props.h"
#include "tetr_elastic_spin_parts.h"
#include "simple_viscoel_spin_parts.h"
#include "simple_viscoel_parts_props.h"
#include "simul_props.h"
#include "elastic_triangle_parts.h"
#include "triangle_parts_props.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

// @Todo : a nicer way to make mesh from config, so that we can keep config from other meshes ?

/*
 Mesh is a wrapper class containing both a (pointer to) wall (*cell_wall) and its (pointer to) properties (*wall_props)
*/

// Defining Mesh class
class Mesh {
    
public:
    // Actual Creator
    Mesh(const YAML::const_iterator );
    Mesh(const YAML::const_iterator , Part_set_props &);
    
    // Dummy creator
    Mesh() {};
    
    // Where the actualy creation occurs
    void Create_mesh(const YAML::const_iterator , Part_set_props &);
    
    // Empty creation of particle set and properties pointers
    Part_set * cell_wall;
    Part_set_props * wall_props;
    
    // Initiation
    void Initiate();
    
};

//#endif
