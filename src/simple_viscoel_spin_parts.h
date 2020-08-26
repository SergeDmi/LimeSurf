/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/

#include "Aboria.h"
using namespace Aboria;
#include "elastic_spin_parts.h"

class Simple_viscoel_set_props;

class Simple_viscoel_part_set : public Elastic_part_set
{
   
public:
    // Dummy creator
    Simple_viscoel_part_set(Simple_viscoel_set_props *);
    
    // Time step for the particle set
    void NextStep(const Simul_props &);
    
    // Compute forces
    void ComputeForces(const Simul_props &);
    
    // Find nearest neighbours
    void GetNeighbours();
    
    // Export the links
    void Export(int );
    
protected :
    // relaxation times of the links
    std::vector<double> link_times;
    
    // Particles set properties
    Simple_viscoel_set_props *prop;
};



