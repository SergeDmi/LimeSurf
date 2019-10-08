/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "tetr_elastic_spin_parts.h"
#include <iostream>
//#include "assert_macro.h"
#include <limits>
#include "tinyply.h"
#include "elastic_parts_props.h"

using namespace tinyply;
using namespace std;

const double PI = boost::math::constants::pi<double>();

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
    std::cout << "# power law : " << power_law << std::endl;
    double k0;
    // Here the status is 0 : not an internal surface
    double status=0.0;
    double dist;
    neigh_pairs pairs_i;
    neigh_pairs pairs_j;
    int si,sj;
    double k_elast=prop->k_elast;
    
    std::cout << "by now we have  " << n_springs << "springs" << std::endl;
    int count=0;
    for(auto const& linker: springs) {
        if (get<4>(linker) ==1.0 ) {
            count++;
        }
    }
    std::cout << "counted linkers in subset " << count << std::endl;
    
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
                
                // Now the easiest part : we create the link
                //double k0=k_elast/(4.0*dist*dist);
                
                switch(power_law) {
                    case 1 : k0=k_elast*dist;           // k_elast is Y :  N/m^2
                        break;
                    case 2 : k0=k_elast/(dist*dist);    // k_elast is Y h : N/m
                        break;
                    case 3 : k0=k_elast/pow(dist,3.0);  // k_elast is Y : N/m^2
                        break;
                }
                // TEMPORARY
                //k0=0.0;
                
                
                if (std::isnan(k0)) {
                    std::cerr << "Error : a link has a NaN k0 " << std::endl;
                } else {
                    n_springs++;
                    link linker{ix,jx,k0,dist,status};
                    springs.push_back(linker);
                }
                
                
                
            }
            
        }
        
    }
    
    std::cout << "by now we have  " << n_springs << "springs" << std::endl;
    
    
    
    
}


// Makes sure everything is in place
void Tetr_elastic_part_set::GetStarted(){
    
    particles.init_id_search();
    Part_set::FindBounds();
    GetNeighbours();
    
    Part_set::CheckPartSet();
    
}

