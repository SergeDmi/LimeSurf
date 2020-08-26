/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "viscous_spin_parts.h"
#include <iostream>
#include <fstream>
#include <limits>
#include "tinyply.h"
#include "viscous_parts_props.h"
using namespace tinyply;
using namespace std;
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

/*
 Old deprecated class : never mind the junk
*/


const double PI = boost::math::constants::pi<double>();

// Dummy constructor
Viscous_part_set::Viscous_part_set() : Part_set() {
    number=0;
    diverging=false;
   
}   

Viscous_part_set::Viscous_part_set(Viscous_set_props * p) : Part_set(p) , prop(p) {
    number=0;
    prop=p;
    double L=prop->L;
    diverging=false;
   
    //double x_conf=prop->x_conf;
    //double y_conf=prop->y_conf;
    //double z_conf=prop->z_conf;
    // Makin da functions yo
    /*
    if (prop->x_conf>0) {
        std::cout << "# There will be confinement along x" << std::endl;
        add_x_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[0]-=pp->x_conf*((posi[0]>pp->x_max)*(posi[0]-pp->x_max) + (posi[0]<-pp->x_max)*(posi[0]+pp->x_max)); };
    } else {
        add_x_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    
    if (prop->y_conf>0) {
        std::cout << "# There will be confinement along y" << std::endl;
        add_y_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[1]-=pp->y_conf*((posi[1]>pp->y_max)*(posi[1]-pp->y_max) + (posi[1]<-pp->y_max)*(posi[1]+pp->y_max)); };
    } else {
        add_y_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    
    if (prop->z_conf>0) {
        std::cout << "# There will be confinement along z" << std::endl;
        add_z_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[2]-=pp->z_conf*((posi[2]>pp->z_max)*(posi[2]-pp->z_max) + (posi[2]<-pp->z_max)*(posi[2]+pp->z_max)); };
    } else {
        add_z_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    */
}

// If particles need to be created from properties
void Viscous_part_set::create() {
    //Part_set::create();
    std::cout << "# starting particle creation ; load " << prop->fname_in << std::endl;
    /*
    if (prop->load_from_file>0) {
        std::string test=prop->fname_in;
        std::cout << "# ... Trying to load " << test << std::endl;
        number=Part_set::load_from_text();
        
    }
    else {
    if (prop->init_shape==0){
        number=PutOnSphere();
    }
    
    if (prop->init_shape==1){
        number=PutOnSheet();
    }
    }
     */
     number=Part_set::load_from_text();
     std::cout << "# created " << number << "particles, with expected R0 " << prop->R0 << std::endl;
     
}



// Loading from ply file
// Copied from tinyply's example
// @TODO : make sure that memcpy is cosher

// Put many particles on a sphere at random
int Viscous_part_set::PutOnSphere(){
    int N=prop->init_number;
    double R=prop->init_radius;
    int i;
    double L=prop->L;
    vdouble3 pos;
    vdouble3 zero(0,0,0);
    vdouble3 cent(L/2.0,L/2.0,L/2.0);
    
    double mindist=prop->minR;
    std::uniform_real_distribution<double> uni(0,1);
    std::default_random_engine generator;

    for ( i = 0; i < N; ++i) {
        int neighbours=1;
        typename particle_type::value_type p;
        double theta;
        double phi;
      
        while (neighbours>0) {
            neighbours=0;
            theta = uni(generator)*2*PI;
            phi = uni(generator)*PI;
            pos=vdouble3(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi));
            get<position>(p) = R*pos+cent;
          
            get<orientation>(p) = pos;
            get<force>(p) = zero;
            get<torque>(p) = zero;
            
            // Using aboria's power
            
            for (auto tpl=euclidean_search(particles.get_query(),get<position>(p),mindist); tpl!=false ; ++tpl) {
                neighbours++;
                break;
            }
             
        }
        get<state>(p) = neighbours;
        particles.push_back(p);
    }
    return i;
}

// Put many particles on a sheet with hexagonal lattice
int Viscous_part_set::PutOnSheet(){
    double R=prop->init_radius;
    double L=prop->L;
    vdouble3 pos1;
    vdouble3 pos2;
    vdouble3 up(0,0,1);
    vdouble3 zero(0,0,0);
    vdouble3 cent(L/2.0-R,L/2.0-R,L/2.0);
    double R0=prop->R0;
    double h=R0*sqrt(3.0)/2.0;
    int Nrows=R/(2*h);
    int Ncols=R/R0;
    int count=0;
    
    typename particle_type::value_type p1;
    typename particle_type::value_type p2;
    for (int i = 0; i < Ncols; ++i) {
        for (int  j = 0; j < Nrows+1; ++j) {
            pos1=vdouble3(i*h*2,j*R0,2*R0*sin((2*PI*i*h*2)/R));
            pos2=vdouble3(i*h*2+h,j*R0+0.5*R0,2*R0*sin((2*PI*(i*h*2+h))/R));
            get<position>(p1) = pos1+cent;
            get<position>(p2) = pos2+cent;
            get<orientation>(p1) = up;
            get<orientation>(p2) = up;
            get<force>(p1) = zero;
            get<torque>(p1) = zero;
            get<force>(p2) = zero;
            get<torque>(p2) = zero;
            count+=2;
            particles.push_back(p1);
            particles.push_back(p2);
        }
    }
    return count;
}
    

// Makes sure everything is in place
void Viscous_part_set::GetStarted(){
    Part_set::GetStarted();
    CheckBoxSize();
    particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
    
}

// Check if box size is OK
// @TODO : repair !
void Viscous_part_set::CheckBoxSize() {
    vdouble3 bottomleft(INFINITY,INFINITY,INFINITY);
    vdouble3 topeuright(-INFINITY,-INFINITY,-INFINITY);
    for (int i = 0; i < number; ++i) {
        vdouble3 posi=get<position>(particles[i]);
        for (int ix=0; ix<3; ++ix) {
            if (posi[ix]<bottomleft[ix]) {
                bottomleft[ix]=posi[ix];
            }
            else if (posi[ix]>topeuright[ix]) {
                topeuright[ix]=posi[ix];
              
            }
        }
    }
    
    for (int ix=0; ix<3; ++ix) {
        double dix=topeuright[ix]-bottomleft[ix];
        bottomleft[ix]=bottomleft[ix]-0.15*dix;
        topeuright[ix]=topeuright[ix]+0.15*dix;
        if (bottomleft[ix]<prop->corner_0[ix]) { prop->corner_0[ix]=bottomleft[ix];}
        if (topeuright[ix]>prop->corner_1[ix]) { prop->corner_1[ix]=topeuright[ix];}
    }
        std::cout << "# bounding box from " << prop->corner_0 << " to " << prop->corner_1 << std::endl;
}


// Get the closest neighbours
void Viscous_part_set::GetNeighbours() {
    int count=0;
    int tnuoc=666;
    int ide;
    int n;
    vdouble3 dir(0,0,0);
    vdouble3 pos1(0,0,0);
    
    for (int i = 0; i < number; ++i) {
        n=0;
        neigh_pairs pairs;
        pos1=get<position>(particles[i]);
        int idi=get<id>(particles[i]);
        //for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
        for (auto tpl = euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax); tpl!=false ; ++tpl ) {
            //const typename particle_type::value_type& j = std::get<0>(tpl);
            //const auto j = std::get<0>(*tpl);
            ide=get<id>(*tpl);
            if (ide!=idi) {
                    dir=get<position>(*tpl)-pos1;
                    pair_n ppp(ide,sqrt(dir.squaredNorm()));
                    pairs.push_back(ppp);
                    n++;
            }
        }

        while (n>prop->max_neighbours) {
            n=pop_furthest_neighbour(&pairs,n);
        }
        get<nn>(particles[i])=n;
        get<neighbours>(particles[i])=pairs;
        
        if (n>count) {
            count=n;
        }
        if (n<tnuoc) {
            tnuoc=n;
        }
            
    }
    
    max_neighbours=count;
  
}

// we should check we got a decent part set after loading from text
void Viscous_part_set::CheckPartSet() {
    int count=0;
    int tnuoc=666;
    vdouble3 orsi;
    vdouble3 zero(0,0,0);
    double normi;
    for (int i = 0; i < number; ++i) {
        orsi=get<orientation>(particles[i]);
        normi=orsi.squaredNorm();
        if (std::isnan(normi)) {
            get<orientation>(particles[i])=zero;
        }
        neigh_pairs pairs=get<neighbours>(particles[i]);
        int n=static_cast<int>(pairs.size());
        
        if (n>count) {
            count=n;
        }
        if (n<tnuoc) {
            tnuoc=n;
        }
    }
    std::cout << "# Maximum neighbours number : " << count << std::endl;
    std::cout << "# Minimum neighbours number : " << tnuoc << std::endl;
}

// removes furthest neighbour from particle
int Viscous_part_set::pop_furthest_neighbour(neigh_pairs * pairs, int n) {
    double dist=0;
    int ex=-1;
    for (int i=0;i<n;++i)
    {
        pair_n ppp=pairs->at(i);
        if (ppp.second>dist) {
            dist=ppp.second;
            ex=i;
        }
    }
    pairs->erase(pairs->begin()+ex);
    return n-1;
}


// One simulation step
void Viscous_part_set::NextStep(const Simul_props & simul_prop){
    ComputeForces(simul_prop);
    AddConfinementForces(simul_prop);
    IntegrateForces(simul_prop);
    //if (!(prop->elastic)) {
    //    particles.update_positions();
    //}
    Part_set::ClearForces();
}


// Computing forces with viscous setting (i.e. changeable nearest neighbours)
//@TODO : verify this part
void Viscous_part_set::ComputeForces(const Simul_props & simul_prop){
    
    vdouble3 posi;
    vdouble3 orsi;
    int idi,idj;
    double p_att=(1.0+prop->p_att)/2.0;
    double p_align=(1.0+prop->p_align)/2.0;
    double p_rep=(prop->p_rep+prop->p_att)/2.0;
    int neibs;
    for (int i = 0; i < number; ++i) {
        neibs=0;
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        idi=get<id>(particles[i]);
        //for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
        for (auto tpl = euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax); tpl!=false ; ++tpl ) {
             //typename particle_type::value_type j = std::get<0>(*tpl);
             //const typename particle_type::value_type& j = std::get<0>(*tpl);
            
            idj=get<id>(*tpl);
            if (idi!=idj) {
                neibs++;
                vdouble3 posj=get<position>(*tpl);
                vdouble3 orsj=get<orientation>(*tpl);
                vdouble3 sumo=orsi+orsj;
                vdouble3 dxij=posj-posi;
                double nsqrij=std::max(dxij.squaredNorm(),prop->minR);
                // Bending torque
                get<torque>(particles[i])-=(prop->k_bend)*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                // Forces : Lennard Jones & alignement
                get<force>(particles[i])+=(-dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att))
                                       +sumo*(0.5*prop->k_align*(dxij.dot(sumo)/pow(nsqrij,p_align))));
            }
        }
        
        get<nn>(particles[i])=neibs;
    }
     //* */
}
