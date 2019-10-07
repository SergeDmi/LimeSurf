#include <random>
#include "Aboria.h"
using namespace Aboria;
//#include "spin_parts_props.h"
#include "simul_props.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>

#ifndef SPIN_PARTS_H
#define SPIN_PARTS_H

class Part_set_props;

class Part_set {
    
    friend class Elastic_part_set;
    friend class Tetr_elastic_part_set;
    friend class Simple_viscoel_part_set;
    friend class Triangle_part_set;
protected:
    typedef std::pair <int,double> pair_n;
    typedef std::vector<pair_n> neigh_pairs;
    typedef std::vector<int> in_faces;
    
    // We decided to use Aboria to create the particle set
    // It should be super convenient if we do a liquid particle set
    // For now, it's way overkill
    
    ABORIA_VARIABLE(orientation,vdouble3,"orientation");
    ABORIA_VARIABLE(torque,vdouble3,"orientation");
    ABORIA_VARIABLE(force ,vdouble3,"orientation");
    ABORIA_VARIABLE(neighbours,neigh_pairs,"neighbours");
    ABORIA_VARIABLE(nn,int,"neighbour number");
    ABORIA_VARIABLE(state,double,"state");
    ABORIA_VARIABLE(nface ,int,"face number");
    ABORIA_VARIABLE(in_triangles,in_faces,"in_triangles");
    
    // The particle set
    //typedef Particles<std::tuple<orientation,neighbours,force,torque,nn,nface,state>,3> particle_type;
    typedef Particles<std::tuple<orientation,neighbours,force,torque,nn,nface,state,in_triangles>,3> particle_type;
    typedef typename particle_type::position position;
    
    // This is convenient if we read faces from a ply file
    typedef struct  {int x,y,z ;} face;
    typedef std::vector<face> face_list;
    
     // This is convenient if we read tetrahedra from a ply file
    typedef struct  {int x,y,z,k;} tetr;
    typedef std::vector<tetr> tetr_list;

public:
    // Creator
    Part_set();
            
    // Creator
    Part_set(Part_set_props*);
    
    // Actually populates the particle set
    void create();
    
    // Properties
    Part_set_props* prop;
    
    // Check divergence of simul
    bool is_diverging() {return diverging ;}
    
    // returns number of particles
    int num() { return number;}
    
    // Performs a time step
    virtual void NextStep(const Simul_props &);
    
    // APply the computed forces
    void IntegrateForces(const Simul_props &);
    
    // Get rids of far away neighbour...
    int pop_furthest_neighbour(neigh_pairs*,int );
    
    // Compute the forces on each particle
    virtual void ComputeForces(const Simul_props &);
    
    // Confining potential
    virtual void AddConfinementForces(const Simul_props &);
    
    // Find the neighbours of every particle
    virtual void GetNeighbours();
    
    // Export to text
    virtual void Export(int);
    
    // Export to ply
    void Export_bly(int,const Simul_props & , double );
    
    // Getting ready to simulate
    virtual void GetStarted();
    
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
    int load_from_text();
    
    // Max number of neighbours
    int max_neighbours;
    
    // Normalizes the normals
    void RenormNorms();
    
    // Find furthest points in any direction
    void FindBounds();
    
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
    double bounds[6];
    
    // Confinement forces
    //std::function<void(vdouble3& ,const vdouble3)> add_x_conf;
    //std::function<void(vdouble3& ,const vdouble3)> add_y_conf;
    //std::function<void(vdouble3& ,const vdouble3)> add_z_conf;
};

#endif
