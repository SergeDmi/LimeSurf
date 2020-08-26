/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include "Aboria.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include <iostream>
#include <fstream>
#include <limits>
#include "tinyply.h"
#include "mesh.h"
using namespace Aboria;


/*
 Mesh is a wrapper class containing both a (pointer to) wall (*cell_wall) and its (pointer to) properties (*wall_props)
*/

Mesh::Mesh(const YAML::const_iterator config )  {
    // wall & props are null pointers to start with
    cell_wall=nullptr;
    wall_props=nullptr;
    // We then create props from the config
    Part_set_props  props =  Part_set_props(config);
    // And then create mesh from the props
    Create_mesh(config,props);
}

// Dummy constructor
Mesh::Mesh(const YAML::const_iterator config , Part_set_props & props) {
    // wall & props are null pointers to start with
    cell_wall=nullptr; 
    wall_props=nullptr;  
    // And then create mesh from the props
    Create_mesh(config,props);
    }
    
// Mesh creation from config and props
void Mesh::Create_mesh(const YAML::const_iterator config , Part_set_props & propal) {
    int mechanics=propal.mechanics;
    // mechanics is the type of mechanics, see Part_set_props
    switch (mechanics) {
        // NOT IMPLEMENTED : viscous simulation
        case 0 : {
            Part_set_props  * props = new Part_set_props(config);
            wall_props=props;
            Part_set * wall=new Part_set(props);
            cell_wall=wall;
            std::cerr << "# NOT IMPLEMENTED : viscous simulation" << std::endl;
            break; }
        // Elastic simulation of thin (2D manifold) mesh
        case 1 : {  
            Elastic_set_props  * props = new Elastic_set_props(config,propal);
            wall_props=props;
            Elastic_part_set * wall=new Elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting 2D elastic simulation" << std::endl;
            break; }
        // NOT IMPLEMENTED Visco-elastic simulation 
        case 2 : {
            Simple_viscoel_set_props *props = new Simple_viscoel_set_props(config);
            wall_props=props;
            Simple_viscoel_part_set * wall = new Simple_viscoel_part_set(props);
            cell_wall=wall;
            std::cerr << "# NOT IMPLEMENTED : visco-elastic simulation" << std::endl;
            break; }
        // NOT IMPLEMENTED : Thin + bending rigidity
        case 3 : {
            Triangle_set_props  * props = new Triangle_set_props(config);
            wall_props=props;
            Triangle_part_set * wall=new Triangle_part_set(props);
            cell_wall=wall;
            std::cerr << "# NOT IMPLEMENTED : 2D elastic simulation + bending" << std::endl;
            break; }                           
        // Thick (3D manifold) elastic shell
        case 4 : {
            Elastic_set_props  * props = new Elastic_set_props(config);
            wall_props=props;
            Tetr_elastic_part_set *wall = new Tetr_elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting 3D elastic simulation" << std::endl;
            break; }
    }
    
}   


// Initiation
void Mesh::Initiate() {
    // Create the cell wall for real
    cell_wall->create();
    // And gets it rolling
    cell_wall->GetStarted();
}

