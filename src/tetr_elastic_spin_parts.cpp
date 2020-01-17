#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "tetr_elastic_spin_parts.h"
#include <iostream>
#include <limits>
#include "tinyply.h"
#include "elastic_parts_props.h"

using namespace tinyply;
using namespace std;

/*
 A class implementing a set of particles connected by springs
Some are on the surface, some are not...
 */

// Dummy constructor
Tetr_elastic_part_set::Tetr_elastic_part_set(Elastic_set_props * p) : Elastic_part_set(p)
{
    
};

// Here we populate the spring set from the face list
void Tetr_elastic_part_set::GetNeighbours() {
    // First we get all the triangles
    Elastic_part_set::GetNeighbours();
    // Now comes the heavy lifting : trying to get the tetrahedrons
    int ix,jx,kx,mx;
    vdouble3 dir,mir,are;
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    int power_law=prop->power_law;
    Part_set::PrintIfVerbose("power law : " + std::to_string(power_law) );
    double k0;
    // Here the status is 0 : not an internal surface
    double status=0.0;
    double dist;
    neigh_pairs pairs_i;
    neigh_pairs pairs_j;
    int si,sj;
    double k_elast=prop->k_elast;
    
    // We go through all the faces and populate interactions
    for (int i = 0; i < n_tetras; ++i) {
        //std::cout << tetrahedra[i].x << " , " << tetrahedra[i].y << " , " << tetrahedra[i].z << " , " << tetrahedra[i].k << " , " << std::endl;
        for (int j=0;j<6;++j) {
            // For all edges if the face
            switch (j) {
                case 0 : ix=tetrahedra[i].x;jx=tetrahedra[i].y ; break ;
                case 1 : ix=tetrahedra[i].y;jx=tetrahedra[i].z ; break ;
                case 2 : ix=tetrahedra[i].z;jx=tetrahedra[i].k ; break ;
                case 3 : ix=tetrahedra[i].k;jx=tetrahedra[i].x ; break ;
                case 4 : ix=tetrahedra[i].x;jx=tetrahedra[i].z ; break ;
                case 5 : ix=tetrahedra[i].y;jx=tetrahedra[i].k ; break ;
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
                dist=sqrt(dir.squaredNorm())/prop->prestrain;
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
                
                k0=k_elast;
                n_springs++;
                //std::cout << "push back linker with ixes " << ix << "," << jx << "    k0=" << k0 << "  dist=" << dist << std::endl;
                link linker{ix,jx,k0,dist,status};
                springs.push_back(linker);
             
                
            }
            
        }
        
    }
    
   Part_set::PrintIfVerbose("by now we have  " + std::to_string(n_springs) + "springs");
       
}


double Tetr_elastic_part_set::ComputeVolumeRatio() {
    
    double tot_volume;
    vdouble3 are;
    vdouble3 posi;
    vdouble3 posj;
    vdouble3 posk;
    vdouble3 posz;
      
      
    // We go through all the faces to compute area
    for (auto const & tetr: tetrahedra) {
        posi=get<position>(particles[tetr.x]);
        posj=get<position>(particles[tetr.y])-posi;
        posk=get<position>(particles[tetr.k])-posi;
        posz=get<position>(particles[tetr.z])-posi;
       
        
        tot_volume+=abs( posj.dot(cross(posk,posz)) );
    }
    tot_volume/=6.0;

    int i,j;
    double sum_l02;
    // We loop over all springs
    // Man C++11 is nice
    // This is the code bottleneck
    for(auto const& linker: springs) {
        sum_l02+=pow(get<3>(linker),2.0);
    }
    
    return tot_volume/sum_l02;
}

void Tetr_elastic_part_set::UpdateLinkerStiffness() {
    
    if (prop->young_modulus>0) {
        double k0=(prop->young_modulus)*(prop->thickness)*(Elastic_part_set::ComputeAreaRatio());

        //std::cout << "computed k_elast = " << k_elast << std::endl;
        switch(prop->power_law) {
            case 1 : k0*=2.0;
                break;
            case 3 : k0*=(6.0/9.0);
                break;
        } 
        for(auto & linker: springs) {
            get<2>(linker)=k0;
        }   
        Part_set::PrintIfVerbose("Changed all links to stiffness " + std::to_string(k0) );
    
    }
    
}


// Makes sure everything is in place
void Tetr_elastic_part_set::GetStarted(){
    Part_set::GetStarted();
    particles.init_id_search();
    Part_set::FindBounds();
    GetNeighbours();
    UpdateLinkerStiffness();
    Part_set::CheckPartSet();
    
}

