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
    double angle,resting_angle;
    int a,b,c,d;
    int ix,iy,iz;
    int j;
    double mean_angle{0};
    double mean_h{0};
    double sum_nrj{0};
    for (auto const & triangle: triangles) {
        a=triangle.x;
        b=triangle.y;
        c=triangle.z;
        
        // just computing h
        vdouble3 ab=get<position>(particles[b])-get<position>(particles[a]);
        vdouble3 ac=get<position>(particles[c])-get<position>(particles[a]);
        
        vdouble3 abXac=cross(ab,ac);
        
        
        for (auto const & other: triangles) {
            
            for ( j=0;j<3;++j) {
                                
                switch (j) {
                    case 0 : ix=other.x ; iy=other.y ; iz=other.z ;  break ;
                    case 1 : ix=other.z ; iy=other.x ; iz=other.y ;  break ;
                    case 2 : ix=other.y ; iy=other.z ; iz=other.x ;  break ;
                }
                if ((ix==b) and (iy==a)) {
               
                    d=iz;
                
                    face_pair faces{a,b,c,d,0.0};
                    angle=ComputeAngle(faces);  
                    mean_angle+=angle;
                    
                    if (prop->imposed_angle>=0) {
                        resting_angle=prop->imposed_angle;
                    } else {
                        resting_angle=angle;
                    }
                    //std::cout << "resting angle :" << resting_angle << std::endl;
                    get<4>(faces)=resting_angle;
                    pairs.push_back(faces);
                    n_pairs++;
                    
                    
                    // Computing total energy
                    
                    
                    vdouble3 ab=get<position>(particles[b])-get<position>(particles[a]);
                    vdouble3 norm1=Part_set::GetNormal(a,b,c);
                    vdouble3 norm2=Part_set::GetNormal(a,d,b);
                    double abn=ab.norm();
                    double n1n=norm1.norm(); // twice the surface area
                    double n2n=norm2.norm(); // twice the surface area
                    double hm=(n1n+n2n)/(2.0*abn);
                    mean_h+=hm; 
                    sum_nrj+=(2.0/3.0)*((n1n+n2n)/4.0)*prop->k_bending*(1.0/2.0)*pow((angle-resting_angle)/hm,2.0);
                    //              surf            curv2 corrections                 curvature squared               
                    
                }
            }
        }
        
    }
    std::cout << "acos(0.5)" << acos(0.5) << "safer : " << safer_acos(0.5) << std::endl;
     std::cout << "mean resting angle :" << mean_angle/n_pairs << std::endl;
     std::cout << "mean h :" << mean_h/n_pairs << std::endl;
     std::cout << "estimated curv : " << pow(prop->curvature_correction,1.0)*mean_angle/mean_h << std::endl;
     std::cout << "estimated radius : " << mean_h/(mean_angle*pow(prop->curvature_correction,1.0)) << std::endl;
     std::cout << "bending nrj : " << pow(prop->curvature_correction,2.0)*sum_nrj << std::endl;
}



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

double Triangle_part_set::ComputeAngle(const vdouble3 & n1, const vdouble3 & n2, const vdouble3 & ab) {    

    return safer_acos(n1.dot(n2)/(n1.norm()*n2.norm()))*sgn(ab.dot(cross(n1,n2)));;
}

double Triangle_part_set::safer_acos(double x) {
  // Safe acos : slow because of branching 
  /*
  if (x < -1.0) x = -1.0 ;
  else if (x > 1.0) x = 1.0 ;
  return acos(x) ;
    */
  
  // Nvidia's : fast (no branching)
  
  double negate = double(x < 0);
  x = abs(x);
  // My own modification !
  // equivalent to min(1.0,x)
  x -= double(x>1.0)*(x-1.0);  
  //x=std::min(x,1.0); <- slower
  //if (x>1.0) {
  //    x=1.0;         <- slower
  //}

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
    ComputeBendingForces(simul_prop);
}

// Bending is super costly !
void Triangle_part_set::ComputeBendingForces(const Simul_props & simul_prop) {
    double angle,resting_angle,abn,n1n,n2n,hm3,pref; //h1,h2,
    int a,b,c,d;
    vdouble3 force1,force2,norm1,norm2,ab,posA,posB;
    // Corrected stiffness
    double stiff=(2.0/3.0)*pow(prop->curvature_correction,2.0)*prop->k_bending;
    // (2.0/3.0) : correction for surface area (one edge has 2/3 the surface area of a face, because triangles)
    // 2.0 : manual correction for curvature. This is more tricky. 
    //std::cout << "-" << stiff << std::flush;
    for (auto const & opair: pairs) {
        a=get<0>(opair);
        b=get<1>(opair);
        c=get<2>(opair);
        d=get<3>(opair);
        resting_angle=get<4>(opair);
        // Slight optimization, significant gains
        posA=get<position>(particles[a]);
        posB=get<position>(particles[b]);
        ab=posB-posA;
        norm1=cross(get<position>(particles[c])-posA,ab);
        norm2=cross(ab,get<position>(particles[d])-posA);
        // Older, cleaner, slower code
        //ab=get<position>(particles[b])-get<position>(particles[a]);
        //norm1=Part_set::GetNormal(a,b,c);
        //norm2=Part_set::GetNormal(a,d,b);
        
        
        abn=ab.norm();
        n1n=norm1.norm(); // twice the surface area
        n2n=norm2.norm(); // twice the surface area
        //angle=ComputeAngle(norm1,norm2,ab);
        // faster version
        angle=safer_acos(norm1.dot(norm2)/(n1n*n2n))*sgn(ab.dot(cross(norm1,norm2)));;
        
        //h1=2.0*n1n/abn;
        //h2=2.0*n2n/abn;
        
        hm3=pow( (n1n+n2n)/(2.0*abn)  ,3.0);
        pref=(stiff*(resting_angle-angle)/(2.0*hm3));
        force1=pref*norm1;
        force2=pref*norm2;
        
        get<force>(particles[c])+=force1;
        get<force>(particles[d])+=force2;
        
        force1=-0.5*(force1+force2);
        get<force>(particles[a])+=force1;
        get<force>(particles[b])+=force1;
        
        /*
        if (angle!=angle) {
            std::cout << "*** angle is nan " << angle << std::endl;
        }
        
        if (norm1.dot(norm2) != norm1.dot(norm2) ) {
            std::cout << "*** n1.n2 nan" << std::endl;
        }
        
        
        if (norm1.dot(norm2) != norm1.dot(norm2) ) {
            std::cout << "*** n1 " << norm1 << "n2 " << norm2 << std::endl;
            std::cout << "*** n1.n2 nan" << std::endl;
        }
    
        if (abn==0.0) {
            std::cout << "zero segment" << std::endl;
        }
        if (norm1.norm()==0.0) {
            std::cout << "zero area1" << std::endl;
        }
        if (norm2.norm()==0.0) {
            std::cout << "zero area2" << std::endl;
        }
        
        if (force1.norm() != force1.norm()) {
            std::cout << "f1" << std::endl;
        }
        
        if (force2.norm() != force2.norm()) {
            std::cout << "f2" << std::endl;
        }
        */
        
    }
}


void Triangle_part_set::GetStarted(){
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

