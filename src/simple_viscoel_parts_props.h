/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <math.h>
#include <sstream>
#include "Aboria.h"
#include "elastic_parts_props.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace Aboria;


#ifndef SIMPLE_VISCOEL_PARTS_PROPS_H
#define SIMPLE_VISCOEL_PARTS_PROPS_H


class Simple_viscoel_set_props : public Elastic_set_props 
{
    friend class Simple_viscoel_part_set;
    
    public:
        /// constructor
        //Simple_viscoel_set_props(const Glossary );
         Simple_viscoel_set_props(const YAML::const_iterator );
        
    protected:
        // distance on which cell wall is plastic
        double size_plastic;
        
        // minimal timescale for plastic effect
        double time_plastic;
};


#endif /* SIMPLE_VISCOEL_PARTS_PROPS_H */

