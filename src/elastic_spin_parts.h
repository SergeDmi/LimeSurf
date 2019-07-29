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

#include "Aboria.h"
using namespace Aboria;
#include "spin_parts.h"


#ifndef ELASTIC_SPIN_PARTS_H
#define ELASTIC_SPIN_PARTS_H

class Elastic_set_props;

class Elastic_part_set : public Part_set
{
    friend class Tetr_elastic_part_set;
    
protected:
    // Link between two points (int,int) with a stiffness and resting length (double double) and status (double)
    typedef std::tuple <int,int,double,double,double> link;
    // Set of all springs
    typedef std::vector<link> spring_set;

public:
    // Dummy creator
    Elastic_part_set(Elastic_set_props *);
    
    // Time step for the particle set
    void NextStep(const Meshless_props*);
    
    // Read from the list of faces to create links and find neighbours
    void GetNeighbours();
    
    // Prepares the simulation
    void GetStarted();
    
    // Compute forces
    void ComputeForces();
    
protected:
    // Springs between particles
    spring_set springs;
  
    // number of springs
    int n_springs;
    
    // Particles set properties
    Elastic_set_props *prop;
    
    // fraction of area per link ; matters for elasticity
    double area_ratio;
    
    // Update area measurements
    void UpdateAreas();
    
    //  mean fraction of area per link  ; matters for pressure
    double mean_area_ratio;
    
    // this computes the module of the force on a linker
    std::function<double(const double,const double,const double)> compute_force;
};
#endif


