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

#ifndef TRIANGLE_PARTS_PROPS_H
#define TRIANGLE_PARTS_PROPS_H


class Triangle_set_props : public Elastic_set_props
{
    friend class Triangle_part_set;

public:
    /// constructor
    Triangle_set_props(const YAML::const_iterator );
    
protected:

    // Bending rigidity (in unit of energy)
    double k_bending;
    
    // A curvature imposed. Warning : complicated for triangulated surface
    // Only validated use is imposed_curvature=0
    double imposed_curvature;
    
    // A parameter to account for the geometry of a triangulated surface
    double curvature_correction;
    
};


#endif /* ELASTIC_PARTS_PROP_H */


