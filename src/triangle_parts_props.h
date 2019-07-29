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
#include "elastic_parts_props.h"

using namespace Aboria;
#ifndef TRIANGLE_PARTS_PROPS_H
#define TRIANGLE_PARTS_PROPS_H


class Triange_set_props : public Elastic_set_props
{
    friend class Triangle_part_set;
    //friend class Tetr_elastic_part_set;
public:
    /// constructor
    Triangle_set_props(const Glossary );
    
    
    
    
protected:
    // prestrain : how much the imported link is stretched
    // double prestrain;
    double rigidity;
};


#endif /* ELASTIC_PARTS_PROP_H */


