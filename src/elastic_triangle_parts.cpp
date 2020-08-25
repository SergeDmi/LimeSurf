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

//#define BENDING_BY_ANGLE 0

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
#ifdef BENDING_BY_ANGLE
    double resting_angle;
#else
    double resting_curv;
#endif
    int v0,v1,vA,vB,vC;
    int found;
    int a,b,c,d;
    int j;
    double mean_angle{0};
    double mean_angle2{0};
    double mean_l{0};
    double mean_hm{0};
    double mean_hm2{0};
    double mean_c{0};
    double mean_c2{0};
    double mean_a{0};
    double sum_nrj{0};
    double sum_nrj_alt{0};
    double mean_radius{0};
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
            mean_angle+=abs(angle);
            mean_angle2+=(angle*angle);
            
            vdouble3 norm1=Part_set::GetNormal(a,b,c);
            vdouble3 norm2=Part_set::GetNormal(a,d,b);
            double abn=ab.norm();
            double n1n=norm1.norm(); // twice the surface area
            double n2n=norm2.norm(); // twice the surface area
            double hm=(n1n+n2n)/(2.0*abn);
            
#ifdef BENDING_BY_ANGLE            
            if (prop->imposed_curvature>=0) {
                resting_angle=prop->imposed_curvature;
            } else {
                resting_angle=angle;
            }
#else 
            if (prop->imposed_curvature>=0) {
                resting_curv=prop->imposed_curvature;
            } else {
                resting_curv=curv_corr*angle/hm;
            }
#endif
            
           
            
            /* Sanity check
            if (norm1.dot(norm2)<0) {
                std::cerr << "OMG opposite triangle orientations" << std::endl;
            }
            */
#ifdef BENDING_BY_ANGLE
            // OLD version : resting angle 
            get<4>(faces)=resting_angle;
#else
            // NEW version : resting curvature
            get<4>(faces)=resting_curv;
#endif
            
            pairs.push_back(faces);
            n_pairs++;
            
            mean_hm+=hm; 
            mean_l+=abn; 
            mean_a+=(n1n+n2n)/(2.0);
            mean_c+=abs(angle)/hm;
            mean_c2+=pow(abs(angle)/hm,2.0);
            mean_hm2+=hm*hm;
            //vdouble3 posmean=(get<position>(particles[a])+get<position>(particles[a])+get<position>(particles[a])/3.0;
            
            
#ifdef BENDING_BY_ANGLE
            // OLD : resting angle 
            sum_nrj+=((n1n+n2n)/2.0)*stiff*(1.0/2.0)*pow(curv_corr*(angle-resting_angle)/hm,2.0);
#else
            // NEW : resting curv
            sum_nrj+=((n1n+n2n)/2.0)*stiff*(1.0/2.0)*pow((curv_corr*angle/hm-resting_curv),2.0);
#endif
            //sum_nrj+=((n1n+n2n)/2.0)*prop->k_bending*(1.0/2.0)*pow((angle/hm-resting_angle/hm),2.0);
            sum_nrj_alt+=prop->k_bending*(1.0/2.0)*pow((angle),2.0);
            //sum_nrj+=(2.0/3.0)*((n1n+n2n)/4.0)*prop->k_bending*(1.0/2.0)*pow((angle-resting_angle)/hm,2.0);
            //              surf            curv2 corrections                 curvature squared               
            
                   
    #ifdef BENDING_BY_ANGLE
            // OLD : resting angle 
            pref=curv_corr*curv_corr*(stiff*(resting_angle-angle)/(2.0*hm*hm*hm));
    #else
            // NEW : resting curv
            pref=curv_corr*(stiff*(resting_curv-curv_corr*angle/hm)/(2.0*hm*hm));
    #endif
            //sum_nrj_alt+=stiff*pow((angle-resting_angle),2.0);

            /*
            vdouble3 force1=pref*norm1;
            vdouble3 force2=pref*norm2;
            
            get<force>(particles[c])+=force1;
            get<force>(particles[d])+=force2;
            
            force1=-0.5*(force1+force2);
            get<force>(particles[a])+=force1;
            get<force>(particles[b])+=force1;
            */    
            
        }    
    }
    /*
    double total_force{0.0};
    vdouble3 curr_force(0,0,0);
    for (auto const part : particles) {
        curr_force=get<force> (part);
        total_force+=curr_force.norm();
    }  
    */
    
    /* Big piece of crap for report.
    std::cout << "mean initial angle :" << mean_angle/n_pairs <<  "mean angle2 : " << mean_angle2/n_pairs << "vs : " << pow(mean_angle/n_pairs,2.0) << std::endl;
    std::cout << "Expected mean angle :" << mean_hm/n_pairs << std::endl;
    std::cout << "mean h :" << mean_hm/n_pairs <<  "   mean h2 :" << mean_hm2/n_pairs << "   vs : " << pow(mean_hm/n_pairs,2.0) << std::endl;
    std::cout << "mean l :" << mean_l/n_pairs << std::endl;
    //std::cout << "estimated curv : " << pow(prop->curvature_correction,1.0)*mean_angle/mean_h << std::endl;
    std::cout << "mean c :" << mean_c/n_pairs  << "estimated curv : " << mean_angle/mean_hm << std::endl;
    //std::cout << "expected corrected energy : " << sum_nrj*pow(2.0/(mean_c/n_pairs),2.0) << std::endl;
    //std::cout << "expected corrected energy2 : " << sum_nrj*(4.0 / (mean_c2/n_pairs) ) << std::endl;
    std::cout << "curv corr. factor : " << std::sqrt((4.0 / (mean_c2/n_pairs) )) << "  - squared : " << (4.0 / (mean_c2/n_pairs) ) << std::endl;
    std::cout << "CORR mean c :" << curv_corr*mean_c/n_pairs  << "CORR mean c2 : " << curv_corr*curv_corr*mean_c2/n_pairs << std::endl;
    
    std::cout << "acos(0.99)=" << acos(0.9999) << "   safer : " << safer_acos(0.9999) << std::endl;
    std::cout << "acos(0.01)=" << acos(0.0001) << "   safer : " << safer_acos(0.0001) << std::endl;
    std::cout << "sum initial angle :" << mean_angle << std::endl;
    std::cout << "n_pairs :" << n_pairs << "   ; n_springs : " << n_springs << std::endl;
    std::cout << "mean initial angle :" << mean_angle/n_pairs <<  "mean angle2 : " << mean_angle2/n_pairs << "vs : " << pow(mean_angle/n_pairs,2.0) << std::endl;
    std::cout << "Expected mean angle :" << mean_hm/n_pairs << std::endl;
    std::cout << "mean h :" << mean_hm/n_pairs <<  "   mean h2 :" << mean_hm2/n_pairs << "   vs : " << pow(mean_hm/n_pairs,2.0) << std::endl;
    std::cout << "mean l :" << mean_l/n_pairs << std::endl;
    //std::cout << "estimated curv : " << pow(prop->curvature_correction,1.0)*mean_angle/mean_h << std::endl;
    std::cout << "mean c :" << mean_c/n_pairs  << "estimated curv : " << mean_angle/mean_hm << std::endl;
    //std::cout << "estimated radius : " << mean_h/(mean_angle*pow(prop->curvature_correction,1.0)) << std::endl;
    std::cout << "estimated radius : " << mean_hm/(mean_angle) << std::endl;
    //std::cout << "bending nrj : " << pow(prop->curvature_correction,2.0)*sum_nrj << std::endl;
    std::cout << "nrj from mean " << prop->k_bending*(1.0/2.0)*mean_a*(mean_angle2/mean_hm2) << "    ;  bending nrj : " << sum_nrj << std::endl;
    std::cout << "mean nrj NOAREA " << prop->k_bending*(1.0/2.0)*n_pairs*(mean_angle2/n_pairs) << "bending nrj NOAREA: " << sum_nrj_alt << std::endl;
    */
     
#ifdef BENDING_BY_ANGLE
    prop->PrintIfVerbose("Bending by angle, not curvature ");
#else
    prop->PrintIfVerbose("Bending by curvature, not angle ");
#endif
    prop->PrintIfVerbose("INITIAL bending nrj : "+std::to_string(sum_nrj));
    //prop->PrintIfVerbose("INITIAL sum of force norm : "+std::to_string(total_force));
    
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
#ifdef BENDING_BY_ANGLE
    // Old
    double resting_angle , hm3 ;
#endif
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
        
#ifdef BENDING_BY_ANGLE
        // OLD : resting angle 
        hm3=pow( hm ,3.0);
        double resting_angle=get<4>(opair);  
        pref=curv_corr*curv_corr*(stiff*(resting_angle-angle)/(2.0*hm3));
#else
        // NEW : resting curvature 
        // resting_curv = get<4>(opair);
        pref=curv_corr*(stiff*(get<4>(opair)-curv_corr*angle/hm)/(2.0*hm*hm));
#endif
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
#ifdef BENDING_BY_ANGLE
    // Old
    double resting_angle , hm3 ;
#else
    // New
    double resting_curv, hm;
#endif
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
        mean_angle+=std::abs(angle);
#ifdef BENDING_BY_ANGLE
        // OLD : resting angle 
        double hm=(n1n+n2n)/(2.0*abn);
        hm3=pow( hm ,3.0);
        resting_angle=get<4>(opair);  
        sum_nrj+=curv_corr_2*((n1n+n2n)/2.0)*pow((angle-resting_angle)/hm,2.0);
#else
        // NEW : resting curv
        hm=(n1n+n2n)/(2.0*abn);
        resting_curv=get<4>(opair);
        sum_nrj+=((n1n+n2n)/2.0)*pow((curv_corr*angle/hm-resting_curv),2.0);
#endif
        //sum_nrj_alt+=pow((angle-resting_angle),2.0);
        sum_nrj_alt+=1.0-norm1.dot(norm2)/(n1n*n2n);

    }
    
    sum_nrj*=0.5*stiff;
    sum_nrj_alt*=0.5*stiff;
    mean_angle/=(1.0*n_pairs);
    prop->PrintIfVerbose("Bending nrj : "+std::to_string(sum_nrj)+"   alt : "+std::to_string(sum_nrj_alt));
    prop->PrintIfVerbose("Mean angle : "+std::to_string(mean_angle));
    return sum_nrj;
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

