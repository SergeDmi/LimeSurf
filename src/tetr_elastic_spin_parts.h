/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   elastic_spin_parts.h
 * Author: dmitrief
 *
 * Created on May 17, 2018, 11:23 AM
 */
//#include <random>
#include "Aboria.h"
using namespace Aboria;
//#include "elastic_parts_props.h"
#include "elastic_spin_parts.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

class Elastic_set_props;

class Tetr_elastic_part_set : public Elastic_part_set
{


public:
    // Dummy creator
    Tetr_elastic_part_set(Elastic_set_props *);
    
    // Read from the list of faces to create links and find neighbours
    void GetNeighbours();
    
    // Prepares the simulation
    void GetStarted();
    
    double ComputeVolumeRatio();
    
    void UpdateLinkerStiffness();
    
};



