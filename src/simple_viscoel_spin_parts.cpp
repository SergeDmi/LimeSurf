/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
//#include "spin_parts.h"
#include "simple_viscoel_spin_parts.h"
#include <iostream>
#include "assert_macro.h"
#include <limits>
#include "simple_viscoel_parts_props.h"



const double PI = boost::math::constants::pi<double>();


// Dummy creator
Simple_viscoel_part_set::Simple_viscoel_part_set(Simple_viscoel_set_props * p) : Elastic_part_set(p), prop(p)
{
    prop=p;
};


// Here we do the actual time step : compute and apply forces
void Simple_viscoel_part_set::NextStep(const Meshless_props* simul_prop){
    // Set all forces & torques to 0
    Part_set::ClearForces();
    // Computing forces and torques
    ComputeForces(simul_prop);
    // Applying the forces
    Part_set::IntegrateForces(simul_prop);
    //From time to time we should check that the normals are normalized
    if (rand()<prop->renorm_rate) {
        RenormNorms();
    }
       
}

void Simple_viscoel_part_set::GetNeighbours() {
    Elastic_part_set::GetNeighbours();
    vdouble3 posij;
    link linker;
    int ni,nj;
    int max_axis;
    
    double dxyz=0;
    double var;
    for (int i=0;i<3;++i) {
        var=bounds[2*i+1]-bounds[2*i];
        if (var>dxyz) {
            dxyz=var;
            max_axis=i;
        }
    }
    
    double ax,px,timescale;
    for (int i=0;i<n_springs;++i){
        linker=springs[i];
        ni=get<0>(linker);   // first vertex
        nj=get<1>(linker);   // second vertex
        //posij=(get<position>(particles[ni])+get<position>(particles[nj]))/2.0;
        ax=(get<position>(particles[ni])[max_axis]+get<position>(particles[nj])[max_axis])/2.0;
        px=std::abs(ax-bounds[2*max_axis+1]);
        ax=std::abs(ax-bounds[2*max_axis]);
        //ax=abs(posij[max_axis]-bounds[2*max_axis]);
        //px=abs(posij[max_axis]-bounds[2*max_axis+1]);
        timescale=prop->time_plastic*(exp(-std::min(ax,px)/prop->size_plastic)) ;
        link_times.push_back(timescale);
    }        
}


// The physics part : computing interaction between vertices
void Simple_viscoel_part_set::ComputeForces(const Meshless_props* simul_prop){
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
    double stretch;
    //double rate=simul_prop->dt/prop->relax;
    double rate;
    // We loop over all springs
    // Man C++11 is nice
    for (int ix=0;ix<n_springs;++ix) {
        link linker=springs[ix];
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
       stretch=(norm2-l0*l0);
       // Force and torque
       get<force>(particles[i])+=k0*dxij*stretch+(orsi*(press*norm2)/ni);
       get<torque>(particles[i])-=align*dir.dot(orsi+orsj)*cross(dir,orsi);
       
       get<force>(particles[j])+=-k0*dxij*stretch+(orsj*(press*norm2)/nj);
       get<torque>(particles[j])-=align*dir.dot(orsi+orsj)*cross(dir,orsj);
       
       get<3>(linker)+=(simul_prop->dt)*link_times[ix]*stretch/l0;
    }
    if (std::isnan(norm2)) {
        // Checking if we diverge
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
}
