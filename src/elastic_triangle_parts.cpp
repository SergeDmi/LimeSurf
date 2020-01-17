/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "elastic_triangle_parts.h"
#include <iostream>
//#include "assert_macro.h"
#include <limits>
#include "tinyply.h"
#include "triangle_parts_props.h"

using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

// Dummy constructor
Triangle_part_set::Triangle_part_set(Triangle_set_props * p) : Elastic_part_set(p), prop(p)
{
    prop=p;
};

// Here we populate the spring set from the face list
void Triangle_part_set::GetNeighbours() {
    Elastic_part_set::GetNeighbours();
    int ix,jx,kx,kkx;
    double spont;
    vdouble3 dir,mir,are;
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    vdouble3 posk(0,0,0);
    double status=1.0;
    double dist;
    neigh_pairs pairs_i;
    neigh_pairs pairs_j;
    int si,sj;
    double k_elast=prop->k_elast;
    double Atot=0;
    double l2tot=0;
    n_springs=0;
    int i,j,l,m;
    int iix;
    //int count;
    in_faces list_faces;
    
    for ( i = 0; i < n_faces; ++i) {
         for ( j=0;j<3;++j) {
             switch (j) {
                 case 0 : ix=triangles[i].x ; break ;
                 case 1 : ix=triangles[i].y ; break ;
                 case 2 : ix=triangles[i].z ; break ;
             }
             
             list_faces=get<in_triangles>(particles[ix]);
             list_faces.push_back(i);
             get<in_triangles>(particles[ix])=list_faces;
             
             
         }
    }
    
    // We go through all the faces and populate interactions
    for ( i = 0; i < n_faces; ++i) {
        for ( j=0;j<3;++j) {
            // For all edges if the face
            switch (j) {
                case 0 : ix=triangles[i].x; jx=triangles[i].y ; kx=triangles[i].z ; break ;
                case 1 : ix=triangles[i].y; jx=triangles[i].z ; kx=triangles[i].x ;break ;
                case 2 : ix=triangles[i].z; jx=triangles[i].x ; kx=triangles[i].y ;break ;
            }
            // Computing triangle surface are
            if (j==0) {
                
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                posk=get<position>(particles[kx]);
                dir=posj-posi;
                mir=posk-posi;
                are=cross(dir,mir);
                Atot+=sqrt(are.squaredNorm());
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
                n_springs++;
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
                
                list_faces=get<in_triangles>(particles[ix]);
                // identifying the other face
                for(auto k : list_faces) {
                    if (k != i) {
                        if (triangles[k].x==ix or triangles[k].y==ix or triangles[k].z==ix) {
                            // We found the other triangle !!!!
                            //l=-1;
                            for ( m=0;j<3;++j) {
                                // For all edges if the face
                                switch (l) {
                                    case 0 : iix=triangles[k].x; break ;
                                    case 1 : iix=triangles[k].y; break ;
                                    case 2 : iix=triangles[k].z; break ;
                                }
                                if (iix != ix and  iix != jx)
                                {
                                    kkx=iix;
                                }
                            }
                            
                        }
                        
                    }
                }
                // Compute the spontaneous curvature !!!!
                
                spont=0.0;
                // Now the easiest part : we create the link
                //double k0=k_elast/(4.0*dist*dist);
                double k0=k_elast/(dist*dist);
                link3 linker{ix,jx,k0,dist,status,kx,kkx,spont};
                l2tot+=dist*dist;
                springs.push_back(linker);
            }
            
        }
        
    }
    Atot/=2.0;
    area_ratio=Atot/l2tot;
    std::cout << "# Found area ratio : " << area_ratio << std::endl;
    std::cout << "# Generated " << n_springs << " springs " << std::endl;
    //particles.erase( std::remove_if(particles.begin(), particles.end(), [](auto& obj){return obj.nfaces == 0;}), particles.end() );
    
    
    /*
     int count=0;
     for (auto it = particles.begin(); it != particles.end(); )
     {
     if (get<nface>(*it)==0) {
     it = particles.erase(it);
     count++;
     }
     else {
     it++;
     }
     }
     number=number-count;
     std::cout << "# Total removed unbound vertices : " << count << std::endl;
     */
    
    
}





// Here we do the actual time step : compute and apply forces
void Triangle_part_set::NextStep(const Simul_props & simul_prop){
    //std::cout << "#" ;
    // Set all forces & torques to 0
    
    // Set all forces & torques to 0
    Part_set::ClearForces();
    
    // Add pressure forces
    // Add pressure forces before confinement to update surface areas !
    Part_set::AddPressureForces(simul_prop);

    // Computing forces and torques
     Elastic_part_set::ComputeForces(simul_prop);
    ComputeForces(simul_prop);
    // Add confinement forces;
    Part_set::AddConfinementForces(simul_prop);
    // Applying the forces
    Part_set::IntegrateForces(simul_prop);

}

// The physics part : computing interaction between vertices
void Triangle_part_set::ComputeForces(const Simul_props & simul_prop){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij,felast;
    double norm2=0;
    vdouble3 dir;
    double p_bend=prop->p_bend/2.0;
    double k_align=prop->k_align*2.0;
    int i,j;
    double l0,k0,project,status;
    double press=simul_prop.pressure*area_ratio;
    
    // We loop over all springs
    // Man C++11 is nice
    for(auto const& linker: springs) {
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        k0=get<2>(linker);  // stiffness
        l0=get<3>(linker);  // resting length
        status=get<4>(linker);  // status
        // We need number of edges to compute pressure
        //ni=get<nn>(particles[i]);
        //nj=get<nn>(particles[j]);
        
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
        //get<force>(particles[i])+=k0*dxij*(norm2-l0*l0)+(orsi*(press*norm2)/ni);
        //get<torque>(particles[i])-=align*dir.dot(orsi+orsj)*cross(dir,orsi);
        
        //get<force>(particles[j])+=-k0*dxij*(norm2-l0*l0)+(orsj*(press*norm2)/nj);
        //get<torque>(particles[j])-=align*dir.dot(orsi+orsj)*cross(dir,orsj);
        felast=k0*dxij*(norm2-l0*l0);
        
        project=k_align*dir.dot(orsi+orsj);
        get<force>(particles[i])+=felast+orsi*(press*norm2)*status;
        get<torque>(particles[i])-=project*cross(dir,orsi)*status;
        
        get<force>(particles[j])+=-felast+orsj*(press*norm2)*status;
        get<torque>(particles[j])-=project*cross(dir,orsj)*status;
        
        
    }
    if (std::isnan(norm2)) {
        // Checking if we diverge
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
}


void Triangle_part_set::GetStarted(){
    Elastic_part_set::GetStarted();
    GetNeighbours();
    
}

