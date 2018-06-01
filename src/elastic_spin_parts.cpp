/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "elastic_spin_parts.h"
#include <iostream>
#include "assert_macro.h"
#include <limits>
//using namespace tinyply;
#include "tinyply.h"
using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

Elastic_part_set::Elastic_part_set(Part_set_props * p) : Part_set(p), prop(p)
//Elastic_part_set::Elastic_part_set(Part_set_props * p)
{
    // Dummy creator
    prop=p;
   
    //particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
}


void Elastic_part_set::NextStep(const Meshless_props* simul_prop){
    Part_set::ComputeForces();
    Part_set::IntegrateForces(simul_prop);
    if (!(prop->elastic)) {
        particles.update_positions();
    }
    Part_set::ClearForces();
}

