/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
#include "glossary.h"
//#include "Aboria.h"
//using namespace Aboria;
#include <sstream>
#include "elastic_parts_props.h"


Elastic_set_props::Elastic_set_props(const Glossary opt) : Part_set_props(opt) {
    
    prestrain=1.0;
    power_law=2;
    opt.set(prestrain, "prestrain");
    opt.set(power_law, "power_law");
    std::cout << "power law : " << power_law << std::endl;
    
}