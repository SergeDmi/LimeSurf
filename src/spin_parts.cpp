#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"

const double PI = boost::math::constants::pi<double>();

Part_set::Part_set(Part_set_props * p) {
    prop=p;
    double L=prop->L;
    //particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false));
    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->minR);
    if (prop->init_shape==0){
        PutOnSphere();
    }

}

int Part_set::num() {return number;}

int Part_set::PutOnSphere(){
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
    //std::random_device r;

    for ( i = 0; i < N; ++i) {
        //bool free_position = false;
        int neighbours=1;
        /*
         * create new particle
         */
        typename particle_type::value_type p;
        double theta;
        double phi;
        
        
        /*
         * randomly choose positions within the domain until one is
         * found with no other particles within a range equal to diameter
         */
        
        
        while (neighbours>0) {
             neighbours=0;
            theta = uni(generator)*2*PI;
            phi = uni(generator)*PI;
            pos=vdouble3(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi));
            get<position>(p) = R*pos+cent;
            //free_position = true;
            get<orientation>(p) = pos;
            get<force>(p) = zero;
            get<torque>(p) = zero;
            /*
             * loop over all neighbouring particles within a euclidean distance
             * of size "diameter"
             */
            for (auto tpl: euclidean_search(particles.get_query(),get<position>(p),mindist)) {
                      neighbours++;
                break;
            }
        }

        get<state>(p) = neighbours;
        particles.push_back(p);
    }
    number=i;
    return i;
}

void Part_set::GetStarted(){
    double L=prop->L;
    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rmax);
    particles.init_id_search();
    GetNeighbours();
}

void Part_set::GetNeighbours() {
    int ide;

    NEIGHBOURS neis;
    int n;
    for (int i = 0; i < number; ++i) {
        n=0;
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
            const typename particle_type::value_type& j = std::get<0>(tpl);
            //int id=get<id>(j);
            //typename particle_type::value_type p=j;
            //std::cout << "the particle id is " << get<id>(p) << std::endl;
            int ide=get<id>(j);
            if (n<MAXNEIGHBOURS) {
                neis[n]=ide;
                n++;
            }
        }
        get<nn>(particles[i])=n;
        get<neighbours>(particles[i])=neis;
    }

}


void Part_set::ViscousStep(const Meshless_props* simul_prop){
    ComputeForcesViscous();
    IntegrateForces(simul_prop);
    ClearForces();
}

void Part_set::ClearForces() {
    vdouble3 zero(0,0,0);
    for (int i = 0; i < number; ++i) {
        get<force>(particles[i])=0;
        get<torque>(particles[i])=0;
    }
    
}

void Part_set::IntegrateForces(const Meshless_props* simul_prop){
    float dt=simul_prop->dt;
    for (int i = 0; i < number; ++i) {
            get<position>(particles[i])+=get<force>(particles[i])*(dt/prop->visco);
            get<orientation>(particles[i])+=cross(get<orientation>(particles[i]),get<torque>(particles[i]))*(dt/prop->Rvisc);
    }
    
}


void Part_set::ComputeForcesViscous(){
    vdouble3 posj;
    vdouble3 orsj;
    vdouble3 posi;
    vdouble3 orsi;
    vdouble3 dxij;
    
    vdouble3 meij;
    double nsqrij;
    int idi;
    int idj;
    double L=prop->L;
    NEIGHBOURS neis;
    double cc_flat;
    double p_att=(1.0+prop->p_att)/2.0;
    double p_rep=(prop->p_rep+prop->p_att)/2.0;
    //particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rmax);
    //particles.init_id_search();
    for (int i = 0; i < number; ++i) {
        idi=get<id>(particles[i]);
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        
        
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
            typename particle_type::value_type j = std::get<0>(tpl);
            idj=get<id>(j);
            
            if (idi<idj) {
                posj=get<position>(j);
                orsj=get<orientation>(j);
                dxij=posj-posi;
                nsqrij=dxij.squaredNorm();
               
                
                
                // Bending torque
                meij=prop->k_bend*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                get<torque>(particles[i])-=meij;
                get<torque>(j)+=meij;

                // Alignement force
                cc_flat=(dxij.dot(orsi+orsj))/(pow(nsqrij,3.5));
                // Lennard-Jones force
                meij=dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att));
                // Adding the forces
                get<force>(particles[i])+=meij+cc_flat*orsi;
                get<force>(j)-=meij+cc_flat*orsj;
                
                
                }

        }
    }
    
}






void Part_set::NextStep(const Meshless_props* simul_prop){
    ViscousStep(simul_prop);
}


void Part_set::Export(int t){
#ifdef HAVE_VTK
    vtkWriteGrid("particles",t,particles.get_grid(true));
#endif
    for (int i = 0; i < number; ++i) {
        //std::cout <<" found "<< i <<std::endl;
        //std::cout <<" Particle  "<< i <<std::endl;
        std::cout <<" "<< get<position>(particles[i])[0] <<" "<< get<position>(particles[i])[1]<<" "<< get<position>(particles[i])[2]<<" "<< get<nn>(particles[i])<<std::endl;
    }
}