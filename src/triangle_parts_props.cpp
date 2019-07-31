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
#include "triangle_parts_props.h"

Triangle_set_props::Triangle_set_props(const Glossary opt) : Elastic_set_props(opt) {
    
    prestrain=1.0;
    
    opt.set(prestrain, "prestrain");
    
    rigidity=0.0;
    
    opt.set(rigidity, "flexural");
}

