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
#include "elastic_parts_props.h"

using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

Elastic_part_set::Elastic_part_set(Elastic_set_props * p) : Part_set(p), prop(p)
{
    prop=p;
    int power_law=prop->power_law;
   
    // We create a function to compute the force to avoid code branches, although it's most likely useless
    switch(power_law) {
   
        case 1 :   compute_force= [] (const double l0, const double norm, const double norm2 ) { return (1.0-l0/norm); };
             break;
         case 2 :  compute_force= [] (const double l0, const double norm, const double norm2 ) { return (norm2-l0*l0); };
             break;
         case 3 :  compute_force= [] (const double l0, const double norm, const double norm2 ) { return (norm2*norm2-norm*l0*l0*l0); };
             break;
     }
};

// Here we populate the spring set from the face list
void Elastic_part_set::GetNeighbours() {
    int ix,jx,kx;
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
            // Computing triangle surface are
            if (j==0) {
                kx=triangles[i].z ;
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                posk=get<position>(particles[kx]);
                dir=posj-posi;
                mir=posk-posi;
                are=cross(dir,mir);
                Atot+=sqrt(are.norm());
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
                switch(power_law) {
                    case 1 : k0=k_elast*dist;           // k_elast is Y :  N/m^2
                        break;
                    case 2 : k0=k_elast/norm2;    // k_elast is Y h : N/m
                        break;
                    case 3 : k0=k_elast/(norm2*dist);  // k_elast is Y : N/m^2 
                        break;
                }
                //k0=k_elast/(dist*dist);
                link linker{ix,jx,k0,dist,status};
                l2tot+=norm2;
                surfaces[ix]+=norm2;
                surfaces[jx]+=norm2;
                springs.push_back(linker);
            }
            
        }
        
    }
    
    
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
    //std::cout << "# Found area ratio : " << area_ratio << std::endl;
    //std::cout << "# Found mean area ratio : " << mean_area_ratio << std::endl;
    //std::cout << "# Found Gamma : " << mean_area_ratio/area_ratio << std::endl;
    //std::cout << "# Generated " << n_springs << " springs " << std::endl;
    //particles.erase( std::remove_if(particles.begin(), particles.end(), [](auto& obj){return obj.nfaces == 0;}), particles.end() );
   
}





// We update the areas and areas ratio !
void Elastic_part_set::UpdateAreas() {
    double Atot=0;
    double l2tot=0;
    double l2mean=0;
    vdouble3 dir,mir,are,dxij;
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    vdouble3 posk(0,0,0);
    std::vector<double> surfaces(number);
    int i,j,ix,jx,kx;
    double dist,norm2,sti;
   
    // We go through all the faces tom compute Atot
    for (int fa = 0; fa < n_faces; ++fa) {
        ix=triangles[fa].x;jx=triangles[fa].y ; kx=triangles[fa].z ;
        posi=get<position>(particles[ix]);
        posj=get<position>(particles[jx]);
        posk=get<position>(particles[kx]);
        dir=posj-posi;
        mir=posk-posi;
        are=cross(dir,mir);
        Atot+=sqrt(are.norm());
        
    }
    
    // We go through all linkers to compute l2tot and surfaces
     for(auto const& linker: springs) {
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        
        // Position & orientation of vertices
        posi=get<position>(particles[i]);
        posj=get<position>(particles[j]);
        
        // Vector & norm of edge
        dxij=posj-posi;
        norm2=dxij.squaredNorm();
        l2tot+=norm2;
        surfaces[ix]+=norm2;
        surfaces[jx]+=norm2;
     }
    
    // we go over all vertices to compute l2mean
     for (int i = 0; i < number; ++i) {
         sti=get<state>(particles[i]);
         if (sti>0.0) {
            l2mean+=surfaces[i]/sti;
         }
     }
     
    Atot/=2.0;
    area_ratio=Atot/l2tot;
    mean_area_ratio=Atot/l2mean;
    //std::cout << "# Found area ratio : " << area_ratio << std::endl;
    //std::cout << "# Found mean area ratio : " << mean_area_ratio << std::endl;
    //particles.erase( std::remove_if(particles.begin(), particles.end(), [](auto& obj){return obj.nfaces == 0;}), particles.end() );
   
}




// Here we do the actual time step : compute and apply forces
void Elastic_part_set::NextStep(const Meshless_props* simul_prop){
    double rand_val;
    //std::cout << "#" ;
    // Set all forces & torques to 0
    Part_set::ClearForces();
    // Computing forces and torques
    ComputeForces();
    // Add confinement forces;
    Part_set::AddConfinementForces();
    // Applying the forces
    Part_set::IntegrateForces(simul_prop);
    //From time to time we should check that the normals are normalized
    rand_val=(double) rand()/RAND_MAX;
    if (rand_val < prop->renorm_rate) {
        RenormNorms();
        UpdateAreas();
    }
       
}

// The physics part : computing interaction between vertices
void Elastic_part_set::ComputeForces(){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij,felast;
    double norm2;
    ///////// WARNING ///////
    // @WARNING
    //  should be area_ratio*2.0 but somehow it is 4.0
    double harat=(area_ratio*4.0);
    // Is it because we're effectively counting twice every linker ?
    // No because that would give just area ratio *1.0
    // OMG I DONT GET IT....

    //double tot_press_force=0;
    vdouble3 dir;
    //double p_bend=prop->p_bend/2.0;
    double k_align=prop->k_align*2.0;
    int i,j;
    double l0,k0,project,status,norm,eff_press,nn_i,nn_j;
    double press=(prop->pressure)*mean_area_ratio;
    //int power_law=prop->power_law;
    //double plaw=1.0*power_law;

    
    // We loop over all springs
    // Man C++11 is nice
    // This is the code bottleneck
    for(auto const& linker: springs) {
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        k0=get<2>(linker);  // stiffness
        l0=get<3>(linker);  // resting length
        status=get<4>(linker);  // status : surface link or not
         
        // Position & orientation of vertices
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        posj=get<position>(particles[j]);
        orsj=get<orientation>(particles[j]);
        
        // Number of neighbours of each vertex (in the inner mesh)
        //nn_i=static_cast<double>(get<state>(particles[i]));
        //nn_j=static_cast<double>(get<state>(particles[j]));
        nn_i=get<state>(particles[i]);
        nn_j=get<state>(particles[j]);
        
        //if (get<state>(particles[i]) != get<nn>(particles[i])) {
        //    std::cerr << "state " << get<state>(particles[i]) <<  "    nn " << get<nn>(particles[i]) << std::endl;
        //    std::cerr << "diffferent state and nn" << std::endl;
        //}
        // Vector & norm of edge
        dxij=posj-posi;
        norm2=dxij.squaredNorm();
        norm=sqrt(norm2);
        dir=dxij/norm;
        
        // Force and torque
        // Force computation is done in a function to avoid branches
        felast=harat*k0*dxij*compute_force(l0,norm,norm2);
        project=k_align*dir.dot(orsi+orsj)*status;
        eff_press=press*norm2*status;
        
        
        
        get<force>(particles[i])+=felast+orsi*eff_press/nn_i;
        get<torque>(particles[i])-=project*cross(dir,orsi);
        
        get<force>(particles[j])+=-felast+orsj*eff_press/nn_j;
        get<torque>(particles[j])-=project*cross(dir,orsj);
        //tot_press_force+=eff_press/nn_i+eff_press/nn_j;
        
    }
    if (std::isnan(norm2)) {
        // Checking if we diverge
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
    //std::cout << " tot press force : " << tot_press_force << std::endl;
}

// Makes sure everything is in place
void Elastic_part_set::GetStarted(){
    //std::cout << "# Getting started in Elastic_part_set " << std::endl;
    //Part_set::GetStarted();
    particles.init_id_search();
    Part_set::FindBounds();
    GetNeighbours();
    
    Part_set::CheckPartSet();
    
}

