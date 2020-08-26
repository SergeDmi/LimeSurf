/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include "Aboria.h"
using namespace Aboria;
#include "elastic_spin_parts.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

/*
 A class implementing a set of particles connected by springs
Some are on the surface, some are not...
 */

class Elastic_set_props;

class Tetr_elastic_part_set : public Elastic_part_set
{


public:
    // Dummy creator
    Tetr_elastic_part_set(Elastic_set_props *);
    
    // Read from the list of faces to create links and find neighbours
    void GetNeighbours();
    
    // Prepares the simulation
    void GetStarted();
    
    double ComputeVolumeRatio();
    
    void UpdateLinkerStiffness();
    
};



