/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
 */
#include <math.h>
#include <sstream>
#include "Aboria.h"
#include "spin_parts_props.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep


using namespace Aboria;
#ifndef ELASTIC_PARTS_PROPS_H
#define ELASTIC_PARTS_PROPS_H


class Elastic_set_props : public Part_set_props 
{
    friend class Elastic_part_set;
    friend class Tetr_elastic_part_set;
    public:
        /// constructor
        //Elastic_set_props(const Glossary );
        Elastic_set_props(const YAML::const_iterator);
        Elastic_set_props(const YAML::const_iterator, Part_set_props & );
        
        void init();
        void Read_config(const YAML::const_iterator);
    
    protected:
        // prestrain : how much the imported link is stretched
        double prestrain;
        int power_law;
        double young_modulus ;      // young modulus of surface
        double thickness ;          // thickness of surface

};


#endif /* ELASTIC_PARTS_PROP_H */

