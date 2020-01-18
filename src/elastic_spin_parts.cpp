#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "elastic_spin_parts.h"
#include <iostream>
#include <limits>
#include "tinyply.h"
#include "elastic_parts_props.h"
using namespace tinyply;

/*
 Elastic_part_set is a set of particles connected by elastic springs
 */
 
// Constructor, from properties p.
Elastic_part_set::Elastic_part_set(Elastic_set_props * p) : Part_set(p), prop(p)
{
    prop=p;
    int power_law=prop->power_law;
   
    // We create a function to compute the force to avoid code branches, although it's most likely useless
    switch(power_law) {
        /*
        case 1 :   compute_force= [] (const double l0, const double norm, const double norm2 ) { return (1.0-l0/norm); };
             break;
         case 2 :  compute_force= [] (const double l0, const double norm, const double norm2 ) { return (norm2-l0*l0); };
             break;
         case 3 :  compute_force= [] (const double l0, const double norm, const double norm2 ) { return (norm2*norm2-norm*l0*l0*l0); };
             break;
        */
        case 1 :  compute_force= [] (const double l0, const double norm2 ) { return (1.0-l0/sqrt(norm2)); };
             break;
        case 2 :  compute_force= [] (const double l0, const double norm2 ) { return (norm2/(l0*l0)-1.0)  ; };
             break;
        case 3 :  compute_force= [] (const double l0, const double norm2 ) { return ( norm2*norm2/(l0*l0*l0*l0)-sqrt(norm2)/l0 ); };
             break;
     }
};

// We compute the area around a linker, defined by the indexes of its two points
double Elastic_part_set::Compute_local_linker_area(int ix, int jx) {
    // The area is the average of the area of the two faces sharing the linker
    double count=0;
    double part_area=0;

    vdouble3 dir,mir,are;
    vdouble3 posi;
    vdouble3 posj;
    vdouble3 posk;

    bool got_points;
   
    // We go through all the faces to compute pressure
    for (auto const & triangle: triangles) {
        got_points=false;
        if (ix==triangle.x or ix==triangle.y or ix==triangle.z) {
            if (jx==triangle.x or jx==triangle.y or jx==triangle.z) {
                got_points=true;
            }
        }
        if (got_points) {
            
            posi=get<position>(particles[triangle.x]);
            posj=get<position>(particles[triangle.y]);
            posk=get<position>(particles[triangle.z]);
            dir=posj-posi;
            mir=posk-posi;
            
            
            are=cross(mir,dir);
            part_area+=are.norm()/2.0;
            count+=1.0;
        }
        
    }
    
    part_area/=count;
    
    return part_area;
    
}

// Here we populate the spring set from the face list
void Elastic_part_set::GetNeighbours() {
    // quite a few variables
    int ix,jx,kx;                   // Indices
    int other_triangle;
    vdouble3 dir,mir,are;           // vectors
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    vdouble3 posk(0,0,0);
    double status=1.0;
    double dist;
    neigh_pairs pairs_i;            // vectors of pairs of neighbours
    neigh_pairs pairs_j;
    int si,sj;
    double k_elast;
    double Atot=0;
    double l2tot=0;
    double l2mean=0;
    double k0,norm2;
    int power_law=prop->power_law;
    n_springs=0;
    std::vector<double> surfaces(number);
    
    // We go through all the faces and populate interactions
    for (int i = 0; i < n_faces; ++i) {
        for (int j=0;j<3;++j) {
            // For all edges if the face
            switch (j) {
                case 0 : ix=triangles[i].x;jx=triangles[i].y ; break ;
                case 1 : ix=triangles[i].y;jx=triangles[i].z ; break ;
                case 2 : ix=triangles[i].z;jx=triangles[i].x ; break ;
            }
            // Computing triangle surface area
            if (j==0) {
                kx=triangles[i].z ;
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                posk=get<position>(particles[kx]);
                dir=posj-posi;
                mir=posk-posi;
                are=cross(dir,mir);
                Atot+=are.norm();
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
                k_elast=prop->k_elast;
                n_springs++;
                // Finding the length of the edge
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                dir=posj-posi;
                // The edge could be pre strained
                norm2=dir.squaredNorm()/pow(prop->prestrain,2.0);
                dist=sqrt(norm2);
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
                if (prop->young_modulus > 0 && prop->thickness > 0) {
                    
                    k_elast=(2.0/3.0)*(prop->young_modulus)*(prop->thickness)*Compute_local_linker_area(ix,jx)/norm2;
                    //std::cout << "computed k_elast = " << k_elast << std::endl;
                    switch(power_law) {
                        case 1 : k0=2.0*k_elast;           
                            break;
                        case 2 : k0=k_elast;        
                            break;
                        case 3 : k0=(6.0/9.0)*k_elast;            
                            break;
                    } 
                }
                else
                {
                    
                    k0=k_elast;
                    /*
                    switch(power_law) {
                        case 1 : k0=k_elast/norm;           // k_elast is Y :  N/m^2
                            break;
                        case 2 : k0=k_elast/norm2;    // k_elast is Y h : N/m
                            break;
                        case 3 : k0=k_elast/(norm2*dist);  // k_elast is Y : N/m^2 
                            break;
                    }
                    */
                }
                
                link linker{ix,jx,k0,dist,status};
                l2tot+=norm2;
                surfaces[ix]+=norm2;
                surfaces[jx]+=norm2;
                springs.push_back(linker);
            }
            
        }
        
    }
    Part_set::PrintIfVerbose("Created  : "+std::to_string(n_springs)+" springs");

    // Here we define particle state as the number of neighbour. 
    // Important to know wether a particle is on the surface or not
     for (int i = 0; i < number; ++i) {
         si=get<nn>(particles[i]);
         if (si>0) {
            l2mean+=surfaces[i]/si;
            get<state>(particles[i])=(double)si;
         }
     }
    Atot/=2.0;
    area_ratio=Atot/l2tot;
    mean_area_ratio=Atot/l2mean;
    Part_set::PrintIfVerbose("Found area ratio : "+std::to_string(area_ratio));
}




// Here we do the actual time step : compute and apply forces
void Elastic_part_set::NextStep(const Simul_props & simul_prop){
    // Set all forces & torques to 0
    Part_set::ClearForces();
    
    // Add pressure forces
    // Add pressure forces before confinement to update surface areas !
    Part_set::AddPressureForces(simul_prop);

    // Computing forces and torques
    ComputeForces(simul_prop);
    
    // Add confinement forces;
    Part_set::AddConfinementForces(simul_prop);
    // Applying the forces
    Part_set::IntegrateForces(simul_prop);
       
}

double Elastic_part_set::ComputeAreaRatio() {
    
    double tot_area;
    vdouble3 are;
    vdouble3 posi;
    vdouble3 posj;
    vdouble3 posk;
      
      
    // We go through all the faces to compute area
    for (auto const & triangle: triangles) {
        posi=get<position>(particles[triangle.x]);
        posj=get<position>(particles[triangle.y]);
        posk=get<position>(particles[triangle.z]);
       
        are=cross(posk-posi,posj-posi);
        tot_area+=are.norm();
    }
    tot_area/=2.0;

    int i,j;
    double sum_l02;
    // We loop over all springs
    // Man C++11 is nice
    // This is the code bottleneck
    for(auto const& linker: springs) {
        if (get<4>(linker)>0) {
            sum_l02+=pow(get<3>(linker),2.0);
        }
    }
    
    return tot_area/sum_l02;
}

// The physics part : computing interaction between vertices
// @TODO : marginal gains in letting compiler to decide memory ?
void Elastic_part_set::ComputeForces(const Simul_props & simul_prop){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij,felast;
    double norm2;

    vdouble3 dir;

    int i,j;
    double k0,l0,norm;
   
    // We loop over all springs
    // Man C++11 is nice
    // This is the code bottleneck
    for(auto const& linker: springs) {
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        k0=get<2>(linker);  // stiffness
        l0=get<3>(linker);  // resting length
        //status=get<4>(linker);  // status : surface link or not
         
        posi=get<position>(particles[i]);
        posj=get<position>(particles[j]);

        // Vector & norm of edge
        dxij=posj-posi;
        norm2=dxij.squaredNorm();
     
     
        // Force computation is done in a function to avoid branches
        felast=k0*dxij*compute_force(l0,norm2); 
        get<force>(particles[i])+=felast;
        get<force>(particles[j])-=felast;
        
    }
     
    if (norm2 != norm2) {
        // Checking if we diverge
        // Faster than diverging=std::isnan(norm2) !
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
    
}

// Makes sure everything is in place
void Elastic_part_set::GetStarted(){
    Part_set::GetStarted();
    particles.init_id_search();
    Part_set::FindBounds();
    GetNeighbours();
    Part_set::CheckPartSet();
    
}

            
