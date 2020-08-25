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
//#include "glossary.h"
#include "Aboria.h"
#include "elastic_parts_props.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace Aboria;
#ifndef TRIANGLE_PARTS_PROPS_H
#define TRIANGLE_PARTS_PROPS_H


class Triangle_set_props : public Elastic_set_props
{
    friend class Triangle_part_set;
    //friend class Tetr_elastic_part_set;
public:
    /// constructor
    //Triangle_set_props(const Glossary );
    Triangle_set_props(const YAML::const_iterator );
    
    
    
protected:
    
    double k_bending;
    
    double imposed_curvature;
    
    double curvature_correction;
    
};


#endif /* ELASTIC_PARTS_PROP_H */


