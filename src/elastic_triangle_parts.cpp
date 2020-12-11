/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include <random>
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "elastic_triangle_parts.h"
#include <iostream>
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

// Here we populate the face pair set from the spring set
// Not performance critical (called in just once per mesh)
void Triangle_part_set::GetNeighbours() {
    double angle;
    double resting_curv;
    int v0,v1,vA,vB,vC;
    int found;
    int a,b,c,d;
    int j;
    double stiff=prop->k_bending/3.0; // (correction for count)
    double curv_corr=prop->curvature_correction;
    double pref;
    
    for(auto const& linker: springs) {
        v0=get<0>(linker);   // first vertex
        v1=get<1>(linker);   // second vertex
        found=0;
        
        for(auto const& triangle: triangles) {
            for ( j=0;j<3;++j) {
                                
                switch (j) {
                    case 0 : vA=triangle.x ; vB=triangle.y ; vC=triangle.z ;  break ;
                    case 1 : vA=triangle.y ; vB=triangle.z ; vC=triangle.x ;  break ;
                    case 2 : vA=triangle.z ; vB=triangle.x ; vC=triangle.y ;  break ;
                }
                
                if ((vA==v0) and (vB==v1)) {
                    a=vA; b=vB; c=vC;
                    found+=1;
                } 
                else if ((vB==v0) and (vA==v1)) {
                    d=vC;
                    found+=1;
                }
                if (found==2) {
                    break;
                }
            }
            if (found==2) {
                break;
            }
        }
        
        if (found==2) {
            // just computing h
            vdouble3 ab=get<position>(particles[b])-get<position>(particles[a]);
            vdouble3 ac=get<position>(particles[c])-get<position>(particles[a]);
            
            face_pair faces{a,b,c,d,0.0};
            angle=ComputeAngle(faces);  
            
            vdouble3 norm1=Part_set::GetNormal(a,b,c);
            vdouble3 norm2=Part_set::GetNormal(a,d,b);
            double abn=ab.norm();
            double n1n=norm1.norm(); // twice the surface area
            double n2n=norm2.norm(); // twice the surface area
            double hm=(n1n+n2n)/(2.0*abn);
            
            if (prop->imposed_curvature>=0) {
                resting_curv=prop->imposed_curvature;
            } else {
                resting_curv=curv_corr*angle/hm;
            }
            
            /* Sanity check
            if (norm1.dot(norm2)<0) {
                std::cerr << "OMG opposite triangle orientations" << std::endl;
            }
            */
            
            get<4>(faces)=resting_curv;
            
            pairs.push_back(faces);
            n_pairs++;
            
        }    
    }
}


// Computes angles between two faces
double Triangle_part_set::ComputeAngle(const face_pair & faces) {
    int a=get<0>(faces);
    int b=get<1>(faces);
    int c=get<2>(faces);
    int d=get<3>(faces);
    vdouble3 ab=get<position>(particles[b])-get<position>(particles[a]);
    vdouble3 ac=get<position>(particles[c])-get<position>(particles[a]);
    vdouble3 da=get<position>(particles[a])-get<position>(particles[d]);
    
    return ComputeAngle( cross(ac,ab) , cross(da,ab) , ab);
}

// Compute angle between two faces defined by their normals, and the edge vector between the two.
// the edge is there to give the sign of the angle, i.e. the faces are ordered.
double Triangle_part_set::ComputeAngle(const vdouble3 & n1, const vdouble3 & n2, const vdouble3 & ab) {    

    return safer_acos(n1.dot(n2)/(n1.norm()*n2.norm()))*sgn(ab.dot(cross(n1,n2)));;
}

// A fast and safe (segfault-free) - yet branchless - way to compute acos.
double Triangle_part_set::safer_acos(double x) const {
  // Safe acos : slow because of branching 
  /*
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acos(x) ;
    */
  
  // Nvidia's : fast (+no branching)
  double negate = double(x < 0);
  x = std::abs(x);
  // My own modification !
  // equivalent to x = min(1.0,x), but no branching !
  x -= double(x>1.0)*(x-1.0);  

  // we let the compiler optimize the algebra
  double ret = -0.0187293;
  ret = ret * x;
  ret = ret + 0.0742610;
  ret = ret * x;
  ret = ret - 0.2121144;
  ret = ret * x;
  ret = ret + 1.5707288;
  ret = ret * std::sqrt(1.0-x);
  ret = ret - 2.0 * negate * ret;
  return negate * 3.14159265358979 + ret;
 
  // Not faster, definitely uglier
  //return negate * 3.14159265358979 + (((((-0.0187293*x)+ 0.0742610)*x - 0.2121144)*x + 1.5707288)* sqrt(1.0-x))*(1.0-2.0*negate);
}

// The physics part : computing interaction between vertices
void Triangle_part_set::ComputeForces(const Simul_props & simul_prop){
    Elastic_part_set::ComputeForces(simul_prop);
    ComputeBendingForces();
}

// Bending is super costly !
void Triangle_part_set::ComputeBendingForces() {
    double angle,abn,n1n,n2n,pref,hm;
    double curv_corr=prop->curvature_correction;
    int a,b,c,d;
    vdouble3 force1,force2,norm1,norm2,ab,posA,posB;
    // Corrected stiffness
    double stiff=prop->k_bending/3.0;
    // This is a correction for geometry : 3 = 2 * 3/2 
    // 2 : because we count 2 triangle areas per edge
    // 3/2 : because there are 3/2 as many edges as faces
    
    // Code bottleneck (for bending)
    for (auto const & opair: pairs) {
        a=get<0>(opair);
        b=get<1>(opair);
        c=get<2>(opair);
        d=get<3>(opair);
       
        // Slight optimization, significant gains by not using get normals but manually computing normals
        posA=get<position>(particles[a]);
        posB=get<position>(particles[b]);
        ab=posB-posA;
        norm1=cross(get<position>(particles[c])-posA,ab);
        norm2=cross(ab,get<position>(particles[d])-posA);
        
        abn=ab.norm();
        n1n=norm1.norm(); // twice the surface area
        n2n=norm2.norm(); // twice the surface area
        
        angle=safer_acos(norm1.dot(norm2)/(n1n*n2n))*sgn(ab.dot(cross(norm1,norm2)));;
        hm=(n1n+n2n)/(2.0*abn);
        
        // NEW : resting curvature
        // resting_curv = get<4>(opair);
        pref=curv_corr*(stiff*(get<4>(opair)-curv_corr*angle/hm)/(2.0*hm*hm));
        force1=pref*norm1;
        force2=pref*norm2;
        
        get<force>(particles[c])+=force1;
        get<force>(particles[d])+=force2;
        
        force1=-0.5*(force1+force2);
        get<force>(particles[a])+=force1;
        get<force>(particles[b])+=force1;
        
    }
    
}

// I don't know wether we actually need this
// But not performance critical
double Triangle_part_set::ComputeBendingEnergy() const {
    double angle,abn,n1n,n2n; //h1,h2,
    double sum_nrj{0},sum_nrj_alt{0};
    // New
    double resting_curv, hm;
    double mean_angle{0};
    int a,b,c,d;
    vdouble3 force1,force2,norm1,norm2,ab,posA,posB;
    double curv_corr=prop->curvature_correction;
    double curv_corr_2=curv_corr*curv_corr;
    double stiff=prop->k_bending/3.0;

    for (auto const & opair: pairs) {
        a=get<0>(opair);
        b=get<1>(opair);
        c=get<2>(opair);
        d=get<3>(opair);
        
       
        // Slight optimization, significant gains
        posA=get<position>(particles[a]);
        posB=get<position>(particles[b]);
        ab=posB-posA;
        norm1=cross(get<position>(particles[c])-posA,ab);
        norm2=cross(ab,get<position>(particles[d])-posA);

        abn=ab.norm();
        n1n=norm1.norm(); // twice the surface area
        n2n=norm2.norm(); // twice the surface area

        angle=safer_acos(norm1.dot(norm2)/(n1n*n2n))*sgn(ab.dot(cross(norm1,norm2)));;
        
        // NEW : resting curv
        hm=(n1n+n2n)/(2.0*abn);
        resting_curv=get<4>(opair);
        sum_nrj+=((n1n+n2n)/2.0)*pow((curv_corr*angle/hm-resting_curv),2.0);

    }
    
    sum_nrj*=0.5*stiff;
    
    return sum_nrj;
}

// Completes properties with derived values if need be
void Triangle_part_set::CompleteProperties() {
    Elastic_part_set::CompleteProperties();
    if (prop->k_bending<0) {
        // Computing rigidity from young modulus
        prop->k_bending=prop->young_modulus*pow(prop->thickness,3.0)/(12.0*(1.0-pow(prop->poisson,2.0)));
        prop->PrintIfVerbose("Warning : Computed k_bending as K = Y h^3 / 12 (1-poisson^2) = "+std::to_string(prop->k_bending));
    }
}

// Getting started before simulation starts
void Triangle_part_set::GetStarted(){
    Triangle_part_set::CompleteProperties();
    Part_set::GetStarted();
    particles.init_id_search();
    Part_set::FindBounds();
    Elastic_part_set::GetNeighbours();
    GetNeighbours();
    Part_set::CheckPartSet();
}


// Here we do the actual time step : compute and apply forces
void Triangle_part_set::NextStep(const Simul_props & simul_prop){
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

