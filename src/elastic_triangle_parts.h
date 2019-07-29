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
#include "spin_parts.h"


#ifndef TRIANGLE_SPIN_PARTS_H
#define TRIANGLE_SPIN_PARTS_H

class Triangle_set_props;

class Triangle_part_set : public Part_set
{
    //friend class Tetr_elastic_part_set;
protected:
    // Link between two points (int,int) with a stiffness and resting length (double double) and status (double)
    // additionally : +2 triangle vertices (int, int) and a spontaneous curvature (double)
    typedef std::tuple <int,int,double,double,double,int,int,double> link3;
    // Set of all springs
    typedef std::vector<link3> spring_set;
    
public:
    // Dummy creator
    Triangle_part_set(Triangle_set_props *);
    
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
    Triangle_set_props *prop;
    
    // fraction of area per link
    double area_ratio;
    
};
#endif
