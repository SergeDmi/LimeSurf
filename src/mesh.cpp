#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include <iostream>
#include <fstream>
//#include "assert_macro.h"
#include <limits>
#include "tinyply.h"
#include "mesh.h"

//using namespace tinyply;
//using namespace std;
//#include "yaml-cpp/yaml.h"  // IWYU pragma: keep


const double PI = boost::math::constants::pi<double>();

Mesh::Mesh(const YAML::const_iterator config )  {
    cell_wall=nullptr;
    wall_props=nullptr;
    Part_set_props  props =  Part_set_props(config);
    
    Create_mesh(config,props);
}

// Dummy constructor
Mesh::Mesh(const YAML::const_iterator config , Part_set_props & props) {
    cell_wall=nullptr; 
    wall_props=nullptr;  
    
    Create_mesh(config,props);
    }
    
void Mesh::Create_mesh(const YAML::const_iterator config , Part_set_props & propal) {
    //propal.Read_config(config);
    int mechanics=propal.mechanics;
    switch (mechanics) {
        case 0 : {
            Part_set_props  * props = new Part_set_props(config);
            wall_props=props;
            Part_set * wall=new Part_set(props);
            cell_wall=wall;
            std::cout << "# Starting viscous simulation" << std::endl;
            break; }
        case 1 : {  
            //Elastic_set_props  * props = new Elastic_set_props(glos);
            Elastic_set_props  * props = new Elastic_set_props(config,propal);
            wall_props=props;
            Elastic_part_set * wall=new Elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting elastic simulation" << std::endl;
            //std::cout << "fname in should in meshless : " << wall_props->fname_in << std::endl;
            break; }
        case 2 : {
            Simple_viscoel_set_props *props = new Simple_viscoel_set_props(config);
            wall_props=props;
            Simple_viscoel_part_set * wall = new Simple_viscoel_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting viscoelastic simulation" << std::endl;
            break; }
        case 3 : {
            Triangle_set_props  * props = new Triangle_set_props(config);
            wall_props=props;
            Triangle_part_set * wall=new Triangle_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting triangular elastic simulation" << std::endl;
            break; }                           
        case 4 : {
            Elastic_set_props  * props = new Elastic_set_props(config);
            wall_props=props;
            Tetr_elastic_part_set *wall = new Tetr_elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting 3D elastic simulation" << std::endl;
            break; }
    }
    
    //node=config;
    //glos=new Glossary;
    std::cout << "# Starting FROM YAML CONFIG" << std::endl;
}   



void Mesh::Initiate() {
    cell_wall->create();
    cell_wall->GetStarted();
}

//void Mesh::Update(YAML::Node config) {
//    wall_props->read_physical_properties(config);
//}
