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
#include "simple_viscoel_parts_props.h"

Simple_viscoel_set_props::Simple_viscoel_set_props(const Glossary opt) : Elastic_set_props(opt) {
    
    time_plastic=DBL_MAX;
    size_plastic=DBL_MAX;
    
    opt.set(time_plastic, "plasticity",0);
    opt.set(size_plastic, "plasticity",1);
    
    
}