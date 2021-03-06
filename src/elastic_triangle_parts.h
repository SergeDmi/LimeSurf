/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include "Aboria.h"
#include "elastic_spin_parts.h"
using namespace Aboria;

#ifndef TRIANGLE_SPIN_PARTS_H
#define TRIANGLE_SPIN_PARTS_H

class Triangle_set_props;

class Triangle_part_set : public Elastic_part_set
{
    //friend class Tetr_elastic_part_set;
protected:
    // A pair of triangles described by 4 point indices (int int int int) as A B C D    
    // And a double (the signed angle) 
    // (With A B the shared edge, C D the opposite vertices)
    typedef std::tuple <int,int,int,int,double> face_pair;
    // Set of all springs
    typedef std::vector<face_pair> pairs_set;
    
public:
    // Dummy creator
    Triangle_part_set(Triangle_set_props *);
    
    // A branchless sign function
    double sgn(double const val) const { return (0.0 < val) - (val < 0.0); }
    
    // A branchless approximated acos function that doesn't crash. Better, faster, stronger.
    double safer_acos(double ) const;
    
    // Time step for the particle set
    void NextStep(const Simul_props &);
    
    // Read from the list of faces to create links and find neighbours
    void GetNeighbours();
    
    // Prepares the simulation
    void GetStarted();
    
    // Compute forces
    void ComputeForces(const Simul_props &);
    
    // Compute bending forces
    void ComputeBendingForces();
    
    // Compute bending energy
    double ComputeBendingEnergy() const;
    
    // Computes angle between two faces
    double ComputeAngle(const face_pair &);
    
    // Computes angle between two faces
    double ComputeAngle(const vdouble3 &, const vdouble3 &, const vdouble3 &);
    
    // Summary of current situation
    virtual void Summary() const { Elastic_part_set::Summary() ; ComputeBendingEnergy(); } ;
        
     // Complete properties
    virtual void CompleteProperties();
    
protected:
    
    // Particles set properties
    Triangle_set_props *prop;
    
    pairs_set pairs;
    
    int n_pairs;
    
    
};
#endif
