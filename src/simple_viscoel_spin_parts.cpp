/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/

#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "simple_viscoel_spin_parts.h"
#include <fstream>
#include <iostream>
#include <limits>
#include "simple_viscoel_parts_props.h"
using namespace std;


const double PI = boost::math::constants::pi<double>();


// Dummy creator
Simple_viscoel_part_set::Simple_viscoel_part_set(Simple_viscoel_set_props * p) : Elastic_part_set(p), prop(p)
{
    std::cerr << "ERROR : Simple visco-elastic set not fully implemented" << std::endl;
    prop=p;
};


// Here we do the actual time step : compute and apply forces
void Simple_viscoel_part_set::NextStep(const Simul_props & simul_prop){
    // Set all forces & torques to 0
    Part_set::ClearForces();
    // Computing forces and torques
    ComputeForces(simul_prop);
    Part_set::AddConfinementForces(simul_prop);
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
        timescale=prop->time_plastic*exp(pow(std::min(ax,px),2.0)/(2.0*pow(prop->size_plastic,2.0)));
        link_times.push_back(timescale);
    }        
}


// The physics part : computing interaction between vertices
void Simple_viscoel_part_set::ComputeForces(const Simul_props & simul_prop){
    vdouble3 posi,posj;
    vdouble3 orsi,orsj;
    vdouble3 dxij;
    int ni,nj;
    double norm2;
    vdouble3 dir;
    vdouble3 tri;
    double align=prop->k_align*2.0;
    int i,j;
    double l0,k0;
    double press=simul_prop.pressure*area_ratio;
    double stretch;

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
       
       get<3>(linker)+=(simul_prop.dt)*stretch/(l0*link_times[ix]);
    }
    if (std::isnan(norm2)) {
        // Checking if we diverge
        std::cerr << "Diverging system " << std::endl;
        diverging=true;
    }
}

// Exporting vertices to a file
void Simple_viscoel_part_set::Export(int t){
    std::string numero(std::to_string(t));
    std::string filename;
    filename="vertices_"+numero+".txt";
    ofstream exportfile;
    exportfile.open(filename);
    //myfile << "Writing this to a file.\n";
    int i,j;
    double l0,k0,t0,norm2,strain;
    vdouble3 posi,posj,dxij;
    exportfile << "# X Y Z X Y Z l0 k0 t0 strain \n";
    for (int ix=0;ix<n_springs;++ix) {
        link linker=springs[ix];
        i=get<0>(linker);   // first vertex
        j=get<1>(linker);   // second vertex
        k0=get<2>(linker);  // stiffness
        l0=get<3>(linker);  // resting length
        t0=link_times[ix];
        posi=get<position>(particles[i]);
        posj=get<position>(particles[j]);
        dxij=posj-posi;
        norm2=dxij.squaredNorm();
        strain=(norm2)/(l0*l0)-1;
        exportfile << posi[0] << " " << posi[1] << " " << posi[2] << " " << posj[0] << " " << posj[1] << " " << posj[2] << " ";
        exportfile << l0 << " " << k0 << " " << t0 << " " << strain << " \n";        
    }
    exportfile.close();
}