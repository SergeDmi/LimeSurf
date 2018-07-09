/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   elastic_parts_prop.h
 * Author: dmitrief
 *
 * Created on July 9, 2018, 11:16 AM
 */
#include <math.h>
#include <sstream>
#include "glossary.h"
#include "Aboria.h"
#include "spin_parts_props.h"

using namespace Aboria;
#ifndef ELASTIC_PARTS_PROP_H
#define ELASTIC_PARTS_PROP_H


class Elastic_set_props : public Part_set_props 
{
    friend class Elastic_part_set;
    
    public:
        /// constructor
        Elastic_set_props(const Glossary );
        

        
        
    protected:
        // prestrain
        double e;
};


#endif /* ELASTIC_PARTS_PROP_H */

