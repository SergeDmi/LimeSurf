#include <random>
#include "Aboria.h"
using namespace Aboria;
//#include "spin_parts_props.h"
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

//#include "glossary.h"


#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

//#ifndef SPIN_PARTS_H
//#define SPIN_PARTS_H

// Todo : a nicer way to make mesh from config, so that we can keep config from other meshes ?

// Mesh is a wrapper class containing both a wall (*cell_wall) and its properties (*wall_props)
class Mesh {
    
    
public:
    // Actual Creator
    Mesh(const YAML::const_iterator );
    Mesh(const YAML::const_iterator , Part_set_props &);
     // Actual Creator
    //Mesh(const Glossary );
    
    Mesh() {};
    
    void Create_mesh(const YAML::const_iterator , Part_set_props &);
    
    //void Read_config(const YAML::const_iterator) {};
    // Empty creation of particle set and properties
    Part_set * cell_wall;
    Part_set_props * wall_props;
    
    // Stores a configuration
    //Glossary glos;
    YAML::const_iterator node;
    
    void Initiate();
    
    //void Update(YAML::Node);
            
};

//#endif
