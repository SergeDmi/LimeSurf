/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>

#include "elastic_spin_parts.h"
#include <iostream>
#include "assert_macro.h"
#include <limits>
#include "tinyply.h"

using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

// Dummy creator
Elastic_part_set::Elastic_part_set(Part_set_props * p) : Part_set(p), prop(p)
{
    prop=p;
}

// Here we populate the spring set from the face list
void Elastic_part_set::GetNeighbours() {
    int ix,jx;
    vdouble3 dir(0,0,0);
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    double dist;
    neigh_pairs pairs_i;
    neigh_pairs pairs_j;
    int si,sj;
    double k_elast=prop->k_elast;
    
    // We go through all the faces and populate interactions
    for (int i = 0; i < n_faces; ++i) {
        for (int j=0;j<3;++j) {
            // For all edges if the face
            switch (j) {
                case 0 : ix=triangles[i].x;jx=triangles[i].y ; break ;
                case 1 : ix=triangles[i].y;jx=triangles[i].z ; break ;
                case 2 : ix=triangles[i].z;jx=triangles[i].x ; break ;
            }
            // Probably useless : noting to how many faces a vertex belong
            get<nface>(particles[ix])+=1;
            // We find if this edge has already been registered !
            bool exists=false;
            pairs_i=get<neighbours>(particles[ix]);
            for(auto const& pair_n: pairs_i) {
                if (pair_n.first==jx) {
                    exists=true;
                }
            }
            
            // If it's a new edge, we register it
            if (exists==false) {
                // Finding the length of the edge
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                dir=posj-posi;
                // The edge could be pre strained
                dist=sqrt(dir.squaredNorm())/prop->e;
                // We create the pair...  Needed for initiation but there could be a better way
                pair_n pair_i(jx,dist);
                pair_n pair_j(ix,dist);
                pairs_j=get<neighbours>(particles[jx]);
                pairs_i.push_back(pair_i);
                pairs_j.push_back(pair_j);
                get<neighbours>(particles[ix])=pairs_i;
                get<neighbours>(particles[jx])=pairs_j;
                
                // This is important ! we count the number of edges for each vertex
                si=static_cast<int>(pairs_i.size());
                sj=static_cast<int>(pairs_j.size());
                get<nn>(particles[ix])=si;
                get<nn>(particles[jx])=sj;
                
                // Now the easiest part : we create the link
                double k0=k_elast/(dist*dist);
                link linker{ix,jx,k0,dist};
                springs.push_back(linker);
            }
            
        }
        
    }
    
}





// Here we do the actual time step : compute and apply forces
void Elastic_part_set::NextStep(const Meshless_props* simul_prop){
    // Set all forces & torques to 0
    Part_set::ClearForces();
    // Computing forces and torques
    ComputeForces();
    // Applying the forces
    Part_set::IntegrateForces(simul_prop);
    //From time to time we should check that the normals are normalized
    if (rand()<prop->renorm_rate) {
        RenormNorms();
    }
       
}

// The physics part : computing interaction between vertices
void Elastic_part_set::ComputeForces(){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij;
    int ni,nj;
    double norm2;
    vdouble3 dir;
    vdouble3 tri;
    double proj;
    double p_bend=prop->p_bend/2.0;
    double align=prop->k_align*2.0;
    int i,j;
    double l0,k0;
    double press=prop->pressure*sqrt(3.0)/4.0;
    
    // We loop over all springs
    // Man C++11 is nice
    for(auto const& linker: springs) {
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        k0=get<2>(linker);  // stiffness
        l0=get<3>(linker);  // resting length
        
        // We need number of edges to compute pressure
        ni=get<nn>(particles[i]);
        nj=get<nn>(particles[j]);
        
        // Position & orientation of vertices
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        posj=get<position>(particles[j]);
        orsj=get<orientation>(particles[j]);

       // Vector & norm of edge
       dxij=posj-posi;
       norm2=dxij.squaredNorm();
       dir=dxij/sqrt(norm2);
       
       // Force and torque
       get<force>(particles[i])+=k0*dxij*(norm2-l0*l0)+(orsi*(press*norm2)/ni);
       get<torque>(particles[i])-=align*dir.dot(orsi+orsj)*cross(dir,orsi);
       
       get<force>(particles[j])+=-k0*dxij*(norm2-l0*l0)+(orsj*(press*norm2)/nj);
       get<torque>(particles[j])-=align*dir.dot(orsi+orsj)*cross(dir,orsj);
    }
    if (std::isnan(norm2)) {
        // Checking if we diverge
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
}

// Makes sure everything is in place
void Elastic_part_set::GetStarted(){
    
    particles.init_id_search();
    GetNeighbours();
    
    Part_set::CheckPartSet();
    
}

