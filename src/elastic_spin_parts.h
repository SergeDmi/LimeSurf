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
#include "spin_parts.h"
//#include "meshless_spin_props.h"
//#include <boost/math/constants/constants.hpp>
//#include <math.h>


class Elastic_part_set : public Part_set
{
    // Link between two points (int,int) with a stifness and resting length (double double)
    typedef std::tuple <int,int,double,double> link;
    typedef std::vector<link> spring_set;
    //typedef std::vector<triang> faces;
 //friend class Part_Set;
    public:
    Elastic_part_set(Part_set_props *);
    Part_set_props *prop; 
    void NextStep(const Meshless_props*);
    void GetNeighbours();
    void GetStarted();
    void ComputeForces();
    spring_set springs;
    
};



