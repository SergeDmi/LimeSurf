/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include "Aboria.h"
using namespace Aboria;
#include "spin_parts.h"
#include "simul_props.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>

#ifndef VISCOUS_SPIN_PARTS_H
#define VISCOUS_SPIN_PARTS_H

/*
 Old deprecated class : never mind the junk
*/

class Viscous_set_props;

class Viscous_part_set : public Part_set {
    
    //friend class Elastic_part_set;
    //friend class Tetr_elastic_part_set;
    //friend class Simple_viscoel_part_set;
    //friend class Triangle_part_set;
protected:
    //typedef std::pair <int,double> pair_n;
    //typedef std::vector<pair_n> neigh_pairs;
    //typedef std::vector<int> in_faces;
    
    // We decided to use Aboria to create the particle set
    // It should be super convenient if we do a liquid particle set
    // For now, it's way overkill
    

public:
    // Creator
    Viscous_part_set();
            
    // Creator
    Viscous_part_set(Viscous_set_props*);
    
    // Actually populates the particle set
    void create();
    
    // Properties
    Viscous_set_props* prop;
    
    
    // Performs a time step
     void NextStep(const Simul_props &);
    
    // Get rids of far away neighbour...
    int pop_furthest_neighbour(neigh_pairs*,int );
    
    // Compute the forces on each particle
    void ComputeForces(const Simul_props &);
    
    // Find the neighbours of every particle
    void GetNeighbours();
    
    
    // Getting ready to simulate
    void GetStarted();
    
    // @TODO : improve CheckBoxSizze()
    // Make sure the box is big enough
    void CheckBoxSize();
    
    // Making sure the Particle set loaded is correct
    void CheckPartSet();
    
    // Set forces & torques to 0
    void ClearForces();
    
    // Put particles on a sphere
    int PutOnSphere();
    
    //Put particles on a sheet
    int PutOnSheet();
    
    // Load particles from file in prop->fname_in
    //int load_from_text();
    
    // Max number of neighbours
    //int max_neighbours;
    
   
    
protected:
    // List of faces
    face_list triangles;
    
    // List of tetrahedra
    tetr_list tetrahedra;
    
    // All particles, aka vertex
    particle_type particles;
    
    // Number of particles
    int number;
    
    // Number of faces
    int n_faces;
    
    // Number of tetrahedra
    int n_tetras;
    
    // Current state of simulation
    bool diverging;

    // Furthest points (xmin xmax ymin ymax zmin zmax)
    //double bounds[6];
    
    // Confinement forces
    //std::function<void(vdouble3& ,const vdouble3)> add_x_conf;
    //std::function<void(vdouble3& ,const vdouble3)> add_y_conf;
    //std::function<void(vdouble3& ,const vdouble3)> add_z_conf;
};

#endif
