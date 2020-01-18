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
    double resting_angle;
    int a,b,c,d;
    int ix,iy,iz;
    int j;
    for (auto const & triangle: triangles) {
        a=triangle.x;
        b=triangle.y;
        c=triangle.z;
        for (auto const & other: triangles) {
            bool got_point=false;
            for ( j=0;j<3;++j) {
                switch (j) {
                    case 0 : ix=other.x ; iy=triangle.y ; iz=triangle.z ;  break ;
                    case 1 : ix=other.z ; iy=triangle.x ; iz=triangle.y ;  break ;
                    case 2 : ix=other.y ; iy=triangle.z ; iz=triangle.x ;  break ;
                }
                if ((ix==b) and (iy==a)) {
                    got_point=true;
                }
                d=iz;
                if (d==c) {
                    std::cerr << "Unexpected behaviour" << std::endl;
                } else {
                    
                    if (got_point) {
                        
                        face_pair faces{a,b,c,d,0.0,0.0};
                        resting_angle=ComputeAngle(faces);
                        get<4>(faces)=resting_angle;
                        pairs.push_back(faces);
                        n_pairs++;
                    }
                    
                }
            }
            
            
           
        }
        
    }
}



double Triangle_part_set::ComputeAngle(const face_pair & faces) {
    int a=get<0>(faces);
    int b=get<1>(faces);
    int c=get<2>(faces);
    int d=get<3>(faces);
    vdouble3 ab=get<position>(particles[b])-get<position>(particles[a]);
    vdouble3 ac=get<position>(particles[c])-get<position>(particles[a]);
    vdouble3 da=get<position>(particles[a])-get<position>(particles[d]);
    
    vdouble3 n1=cross(ac,ab);
    vdouble3 n2=cross(da,ab);
    vdouble3 crossed=cross(n1,n2);
    
    if (n1.dot(n2)<0) {
        std::cout << "weird :::: " << std::endl;
    }
    
    return ComputeAngle(n1,n2,ab);

}

double Triangle_part_set::ComputeAngle(const vdouble3 & n1, const vdouble3 & n2, const vdouble3 & ab) {
    return acos(n1.dot(n2))*sgn(ab.dot(cross(n1,n2)))/(n1.norm()*n2.norm());
}

// The physics part : computing interaction between vertices
void Triangle_part_set::ComputeForces(const Simul_props & simul_prop){
    Elastic_part_set::ComputeForces(simul_prop);
    ComputeBendingForces(simul_prop);
}

void Triangle_part_set::ComputeBendingForces(const Simul_props & simul_prop) {
    double angle,resting_angle;
    int a,b,c,d,h1,h2,abn;
    vdouble3 force1,force2,norm1,norm2,ab;
    double stiff=prop->k_bending;
    for (auto const & opair: pairs) {
    //for (auto const & opair: pairs_set) {
        //angle=ComputeAngle(faces);
        a=get<0>(opair);
        b=get<1>(opair);
        c=get<2>(opair);
        d=get<3>(opair);
        resting_angle=get<4>(opair);
        ab=get<position>(particles[b])-get<position>(particles[a]);
        norm1=Part_set::GetNormal(a,b,c);
        norm2=Part_set::GetNormal(a,c,b);
        angle=ComputeAngle(norm1,norm2,ab);
        abn=ab.norm();
        h1=2.0*norm1.norm()/abn;
        h2=2.0*norm2.norm()/abn;
        force1=(stiff*(resting_angle-angle)/(2.0*h1*h1*h1))*norm1;
        force2=(stiff*(resting_angle-angle)/(2.0*h2*h2*h2))*norm2;
        
        get<force>(particles[c])+=force1;
        get<force>(particles[d])+=force2;
        get<force>(particles[a])-=2.0*(force1+force2);
        get<force>(particles[b])-=2.0*(force1+force2);
    }
}


void Triangle_part_set::GetStarted(){
    Elastic_part_set::GetStarted();
    GetNeighbours();
    
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

