#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include "assert_macro.h"
const double PI = boost::math::constants::pi<double>();

Part_set::Part_set(Part_set_props * p) {
    prop=p;
    double L=prop->L;
    //particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false));
    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rsearch);
    if (prop->init_shape==0){
        number=PutOnSphere();
    }
    
    if (prop->init_shape==1){
        number=PutOnSheet();
    }

}

double Part_set::Reff2(double R) {
    return std::max(R,prop->minR2);
}

double Part_set::Reff(double R) {
    return std::max(R,prop->minR);
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
        typename particle_type::value_type p;
        double theta;
        double phi;
      
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
    
    return i;
}


int Part_set::PutOnSheet(){
    int N=prop->init_number;
    double R=prop->init_radius;
    
    double L=prop->L;
    vdouble3 pos1;
    vdouble3 pos2;
    vdouble3 up(1,0,0);
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
            pos1=vdouble3(i*h*2,j*R0,0);
            pos2=vdouble3(i*h*2+h,j*R0+0.5*R0,0);
            get<position>(p1) = pos1+cent;
            get<position>(p2) = pos2+cent;
            //free_position = true;
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
    
    std::cout << "# created " << count << "particles" << "expected R0 " << prop->R0 << std::endl;
    number=count;
    return count;
}


void Part_set::GetStarted(){
    double L=prop->L;
    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rsearch);
    std::cout << "# initiated neighbour serch xith Rmax" << prop->Rmax << std::endl;
    particles.init_id_search();
    //GetNeighbours();
}

void Part_set::GetNeighbours() {
    int ide;

    NEIGHBOURS neis;
    int n;
    for (int i = 0; i < number; ++i) {
        n=0;
        int idi=get<id>(particles[i]);
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
        //typename particle_type::value_type ii=*particles.get_query().find(idi);
        
        get<nn>(particles[i])=n;
        get<neighbours>(particles[i])=neis;
    }

}


void Part_set::ViscousStep(const Meshless_props* simul_prop){
    
    ComputeForcesViscousFast();
    //GetNeighbours();
    //IntegrateForces(simul_prop);
    IntegrateForcesFast(simul_prop);
    particles.update_positions();
    ClearForces();
}

void Part_set::ClearForces() {
    vdouble3 zero(0,0,0);
    for (int i = 0; i < number; ++i) {
        get<force>(particles[i])=zero;
        get<torque>(particles[i])=zero;
    }
    
}

void Part_set::IntegrateForces(const Meshless_props* simul_prop){
    float dt=simul_prop->dt;
    for (int i = 0; i < number; ++i) {
        //if (i==5) {
        //    std::cout << "# position 5 was " <<  get<position>(particles[i]) << std::endl;
        //    std::cout << "# force on 5 is " <<  get<force>(particles[i]) << std::endl;
        //    std::cout << "# da foooorc is " <<  (get<force>(particles[i]))*(dt/prop->visco) << std::endl;
        //
        //}
            get<position>(particles[i])+=(get<force>(particles[i]))*(dt/prop->visco);
        //if (i==5) {
        //    std::cout << "# position 5 is  " <<  get<position>(particles[i]) << std::endl;
        //}
      
            get<orientation>(particles[i])+=cross(get<orientation>(particles[i]),get<torque>(particles[i]))*(dt/prop->Rvisc);
    }
    
}



void Part_set::IntegrateForcesFast(const Meshless_props* simul_prop){
    float dt=simul_prop->dt;
    Label<0, particle_type> a(particles);
    double Fmax=prop->Fmax;
    Symbol<position> p;
    Symbol<force> f;
    /*
    if (norm(f[a])>Fmax) {
        f[a]*=Fmax/norm(f[a]);
    }
     */
    p[a]+=(dt/prop->visco)*(f[a]);
    
}


void Part_set::ComputeForcesViscous(){
    //vdouble3 posj;
   // vdouble3 orsj;
    vdouble3 posi;
    vdouble3 orsi;
    //vdouble3 dxij;
    
    //vdouble3 meij;
    //double nsqrij;
    int idi;
    //int idj;
    double L=prop->L;
    NEIGHBOURS neis;
    double cc_flat;
    //int n=0;
    double p_att=(1.0+prop->p_att)/2.0;
    double p_rep=(prop->p_rep+prop->p_att)/2.0;
    for (int i = 0; i < number; ++i) {
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
             typename particle_type::value_type j = std::get<0>(tpl);
            vdouble3 posj=get<position>(j);
            vdouble3 orsj=get<orientation>(j);
            vdouble3 sumo=orsi+orsj;
            vdouble3 dxij=posj-posi;
            double nsqrij=std::max(dxij.squaredNorm(),prop->minR);
            // Bending torque
            get<torque>(particles[i])-=(prop->k_bend)*cross(orsi,orsj)/(pow(nsqrij,3.0));;
            // Alignement & lennard Jones force
            //get<force>(particles[i])+=((0.5*prop->k_align)*((dxij.dot(sumo))/(pow(nsqrij,3.5)))*sumo-dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att)));
            get<force>(particles[i])+=(-dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att)));
            
        }
    }
}

void Part_set::ComputeForcesViscousFast(){
    //vdouble3 posj;
    // vdouble3 orsj;
    vdouble3 posi;
    vdouble3 orsi;
    //vdouble3 dxij;
    
    //vdouble3 meij;
    //double nsqrij;
    int idi;
    //int idj;
    double L=prop->L;
    NEIGHBOURS neis;
    double cc_flat;
    //int n=0;
    double p_att=(1.0+prop->p_att);
    double p_rep=(1.0+prop->p_rep);
    Symbol<position> p;
    Symbol<orientation> v;
    Symbol<force> f;
    Symbol<id> id_;
    Label<0, particle_type> a(particles);
    Label<1, particle_type> b(particles);
    auto dx = create_dx(a, b);
    AccumulateWithinDistance<std::plus<vdouble3>> sum(prop->Rmax);
    
    f[a] += (
                  // spring force between particles
                                sum(b, if_else(id_[a] != id_[b],
                                               (-dx*((prop->k_rep)/pow(norm(dx),p_rep)
                                                     -(prop->k_att)/pow(norm(dx),p_att))),0.0)
                                               )
             );
                  //                             -k * (diameter / norm(dx) - 1), 0.0) *
                  //                  dx) /
                  //              mass);
    
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
        //std::cout <<" "<< get<position>(particles[i]) <<" " << get<force>(particles[i]) <<" "<< get<nn>(particles[i])<<" " << get<id>(particles[i])<<std::endl;
        
        std::cout <<" "<< get<position>(particles[i])[0] <<" "<< get<position>(particles[i])[1] << " "<< get<position>(particles[i])[2] << " " << get<force>(particles[i])[0] << " " << get<force>(particles[i])[1] << " " << get<force>(particles[i])[2] <<" "<< get<nn>(particles[i])<<" " << get<id>(particles[i])<<std::endl;
        
    }
}


    //particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rmax);
    //particles.init_id _search();
    //Symbol<position> p;
    //Symbol<velocity> v;
    //Symbol<id> id_;
    //Label<0, container_type> a(particles);
    //Label<1, container_type> b(particles);
    
    //auto dx = create_dx(a, b);
    //AccumulateWithinDistance<std::plus<vdouble3>> sum(diameter);
    //AccumulateWithinDistance<std::plus<vdouble2>> sum(diameter);
    
    //v[a] += dt * (
                  // spring force between particles
    //              sum(b, if_else(id_[a] != id_[b],
    //                             -k * (diameter / norm(dx) - 1), 0.0) *
    //                  dx) /
    //              mass);
    //p[a] += dt * v[a];
    
    /*
    for (int i = 0; i < number; ++i) {
        idi=get<id>(particles[i]);
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        
        
        
        //std::cout <<"# Particle "<< i << " with ID" << idi << std::endl;
        //n=0;
        
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
            typename particle_type::value_type j = std::get<0>(tpl);
            //const typename particle_type::value_type& jj = std::get<0>(tpl);
            // By definition here (?), jj should be const. This is we need to get ugly
            idj=get<id>(j);
            //typename particle_type::value_type& j = jj;
            
            
            //n++;
            if (idj<idi) {
                //j=getParticleByID(idj);
                 //typename particle_type::value_type& j=particles.get_query().find(idj);
                std::cout <<"#           interacts with ID" << idj << std::endl;
                posj=get<position>(j);
                orsj=get<orientation>(j);
                vdouble3 sumo=orsi+orsj;
                dxij=posj-posi;
                //std::cout <<"# dxij " << dxij << std::endl;
                nsqrij=std::max(dxij.squaredNorm(),prop->minR);

                // Bending torque
                meij=prop->k_bend*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                get<torque>(particles[i])-=meij;
                get<torque>(j)+=meij;

                // Alignement force
                cc_flat=(dxij.dot(orsi+orsj))/(pow(nsqrij,3.5));
                // Lennard-Jones force
                meij=dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att));
                std::cout <<"# force " << meij << std::endl;
                //std::cout <<"# else " << cc_flat*orsi << std::endl;
                // Adding the forces
                // For some reasons, much bigger forces on y for a regular lattice... WTF ?
                std::cout <<"# force on " << idi << " was" <<  get<force>(particles[i]) << std::endl;
                //get<force>(particles[i])+=meij+cc_flat*orsi;
                get<force>(particles[i])+=meij;
                std::cout <<"# force on " << idi << "  is" <<  get<force>(particles[i]) << std::endl;
                std::cout <<"# *force on " << idj << " was" << get<force>(j) << std::endl;
                //get<force>(j)-=meij+cc_flat*orsj;
                get<force>(j)-=meij;
                std::cout <<"# *force on " << idj << "  is" << get<force>(j) << std::endl;
                
                //if (get<id>(particles[i])~=idi){
                //    std::cout <<"# counted " << n << " neighbours" << std::endl;
                //}
                assert_true(get<id>(particles[i])==idi);
                assert_true(get<id>(j)==idj);
            }

        }
     
        //int neis=   get<nn>(particles[i]);
        //std::cout <<"# counted " << n << " neighbours" << std::endl;
        //std::cout <<"# was expecting " << neis << " neighbours" << std::endl;
    }
    
     */






