/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   simple_viscoel_parts_props.h
 * Author: dmitrief
 *
 * Created on July 9, 2018, 3:01 PM
 */
#include <math.h>
#include <sstream>
#include "glossary.h"
#include "Aboria.h"
#include "elastic_parts_props.h"

using namespace Aboria;


#ifndef SIMPLE_VISCOEL_PARTS_PROPS_H
#define SIMPLE_VISCOEL_PARTS_PROPS_H


class Simple_viscoel_set_props : public Elastic_set_props 
{
    friend class Simple_viscoel_part_set;
    
    public:
        /// constructor
        Simple_viscoel_set_props(const Glossary );
        
    protected:
        // distance on which cell wall is plastic
        double size_plastic;
        
        // minimal timescale for plastic effect
        double time_plastic;
};


#endif /* SIMPLE_VISCOEL_PARTS_PROPS_H */

