/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include "Aboria.h"
using namespace Aboria;
#include "spin_parts.h"

#ifndef ELASTIC_SPIN_PARTS_H
#define ELASTIC_SPIN_PARTS_H

class Elastic_set_props;

class Elastic_part_set : public Part_set
{
    friend class Tetr_elastic_part_set;
    friend class Triangle_part_set;
    
protected:
    // Link between two points (int,int) with a stiffness and resting length (double double) and status (double)
    typedef std::tuple <int,int,double,double,double> link;
    // Set of all springs
    typedef std::vector<link> spring_set;

public:
    // Dummy creator
    Elastic_part_set(Elastic_set_props *);
    
    // Time step for the particle set
    virtual  void NextStep(const Simul_props &);
    
    // Read from the list of faces to create links and find neighbours
    virtual void GetNeighbours();
    
    // Prepares the simulation
    virtual  void GetStarted();
    
    // Compute forces
    virtual  void ComputeForces(const Simul_props &);
    
    // Compute area associated to two points (i.e. linker extremities)
    double Compute_local_linker_area(int, int);
    
    // Computes something
    double ComputeAreaRatio();
    
    // Computes stretching energy
    double ComputeStretchingEnergy() const;
    
    // Summary of current situation
    virtual void Summary() const { Part_set::Summary() ;  ComputeStretchingEnergy(); } ;
    
protected:
    // Springs between particles
    spring_set springs;
  
    // number of springs
    int n_springs;
    
    // Particles set properties
    Elastic_set_props *prop;
    
    // fraction of area per link ; matters for elasticity
    double area_ratio;
    
    //  mean fraction of area per link  ; matters for pressure
    double mean_area_ratio;
    
    // this computes the module of the force on a linker
    std::function<double(const double,const double)> compute_force;
    
    // this computes the deformation energy of a linker
    std::function<double(const double,const double)> compute_energy;
    
    // Number of pairs of faces to compute angles for rigidity
    int n_pairs;
    
};
#endif


