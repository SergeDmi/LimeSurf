/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <random>
//#include "Aboria.h"
//using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>

#include "elastic_spin_parts.h"
#include <iostream>
#include "assert_macro.h"
#include <limits>
//using namespace tinyply;
#include "tinyply.h"
//#include "spin_parts.h"
using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

Elastic_part_set::Elastic_part_set(Part_set_props * p) : Part_set(p), prop(p)
//Elastic_part_set::Elastic_part_set(Part_set_props * p)
{
    // Dummy creator
    prop=p;
   
    //particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
}


void Elastic_part_set::NextStep(const Meshless_props* simul_prop){
    Part_set::ClearForces();
    ComputeForces();
    Part_set::IntegrateForces(simul_prop);
   
    
    if (rand()<prop->renorm_rate) {
        RenormNorms();
    }
       
}

void Elastic_part_set::ComputeForces(){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij;
    int ni,nj;
    double k_elast=prop->k_elast;
    double norm2;
    vdouble3 dir;
    vdouble3 tri;
    double proj;
    double p_bend=prop->p_bend/2.0;
    int count;
    double align=prop->k_align*2.0;
    int i,j;
    double l0;
    //for (int i = 0; i < number; ++i) {
        //std::cout << "# elastic pot " << k_elast << std::endl;
    for(auto const& linker: springs) {
        i=get<0>(linker);
        j=get<1>(linker);
        l0=get<2>(linker);

        ni=get<nn>(particles[i]);
        nj=get<nn>(particles[j]);
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
         posj=get<position>(particles[j]);
        orsj=get<orientation>(particles[j]);

        dxij=posj-posi;
       norm2=dxij.squaredNorm();
       // We need to make it superlinear
       dir=dxij/sqrt(norm2);
       //get<force>(particles[i])+=k_elast*(dxij-l0*dir);
       get<force>(particles[i])+=k_elast*dir*(norm2-l0*l0)+(orsi*(prop->pressure*norm2)/ni);
       get<torque>(particles[i])-=align*dir.dot(orsi+orsj)*cross(dir,orsi);
       get<force>(particles[j])+=-k_elast*dir*(norm2-l0*l0)+(orsj*(prop->pressure*norm2)/nj);
       get<torque>(particles[j])-=align*dir.dot(orsi+orsj)*cross(dir,orsj);
       if (std::isnan(norm2)) {
           std::cerr << "Diverging system " << std::endl;
           diverging=true;
       }
    }
    //std::cout << "# computed " << count << " neighbours" <<std::endl;
}

// Makes sure everything is in place
void Elastic_part_set::GetStarted(){
    //double L=prop->L;
    //Part_set::CheckBoxSize();
    //particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    //std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
    //particles.update_positions();
    if (prop->elastic) {
        particles.init_id_search();
        GetNeighbours();
    }
    Part_set::CheckPartSet();
    //particles.init_id_search();
}


// Get the closest neighbours
void Elastic_part_set::GetNeighbours() {
    // Reads directly from the ply file
    // not necessarily super fast, but needs to get done only once
    int ix,jx;
    vdouble3 dir(0,0,0);
    vdouble3 posi(0,0,0);
    vdouble3 posj(0,0,0);
    double dist;
    neigh_pairs pairs_i;
    neigh_pairs pairs_j;
    int count=0;
    int tnuoc=666;
    int si,sj;
    int min_id;
    int max_id;
    
    
    for (int i = 0; i < n_faces; ++i) {
        //std::cout<< " new face with " << triangles[i].x << " " << triangles[i].y << " " << triangles[i].z << std::endl;
        for (int j=0;j<3;++j) {
            switch (j) {
                case 0 : ix=triangles[i].x;jx=triangles[i].y ; break ;
                case 1 : ix=triangles[i].y;jx=triangles[i].z ; break ;
                case 2 : ix=triangles[i].z;jx=triangles[i].x ; break ;
                }
            get<nface>(particles[ix])+=1;
            bool exists=false;
            pairs_i=get<neighbours>(particles[ix]);
            //std::cout << "looking for " <<jx << " in " <<std::flush;
            for(auto const& pair_n: pairs_i) {
                //std::cout << pair_n.first << " ";
                if (pair_n.first==jx) {
                   exists=true;
                }
            }
            //std::cout << "result " << exists << std::flush;
               
            if (exists==false) {
                //std::cout << " result " << exists << std::flush;
                //std::cout << " ix  = " << ix << " & jx = " << jx << std::endl;
                posi=get<position>(particles[ix]);
                posj=get<position>(particles[jx]);
                dir=posj-posi;
                dist=sqrt(dir.squaredNorm())/prop->e;
                pair_n pair_i(jx,dist);
                pair_n pair_j(ix,dist);
                
                
                pairs_j=get<neighbours>(particles[jx]);
                pairs_i.push_back(pair_i);
                pairs_j.push_back(pair_j);
                get<neighbours>(particles[ix])=pairs_i;
                get<neighbours>(particles[jx])=pairs_j;
                si=static_cast<int>(pairs_i.size());
                sj=static_cast<int>(pairs_j.size());
                get<nn>(particles[ix])=si;
                get<nn>(particles[jx])=sj;
                link linker{ix,jx,dist};
                //get<0>(linker)=ix;get<1>(linker)=jx;get<2>(linker)=dist;
                springs.push_back(linker);
            }
            //std::cout << std::endl;
        }
        
    }
  
}