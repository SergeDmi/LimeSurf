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
        PutOnSphere(prop->init_number,1.0);
    }

}

int Part_set::num() {return number;}

int Part_set::PutOnSphere(int N, double R){
    int i;
    vdouble3 pos;
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
        vdouble3 zero(0,0,0);
        
        /*
         * randomly choose positions within the domain until one is
         * found with no other particles within a range equal to diameter
         */
        
        
        while (neighbours>0) {
             neighbours=0;
            theta = uni(generator)*2*PI;
            phi = uni(generator)*PI;
            pos=vdouble3(cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi));
            get<position>(p) = pos+2.0*vdouble3(1.0,1.0,1.0);
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
    int ide;
    double L=prop->L;
    NEIGHBOURS neis;
    int n;
    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rmax);
    particles.init_id_search();
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
    float dt=simul_prop->dt;
    /*
     * create symbols and labels in order to use the expression API
     Symbol<position> p;
     Symbol<orientation> v;
     Symbol<id> id_;
     Label<0, particle_type> a(particles);
     Label<1, particle_type> b(particles);
     */
    
    /*
     * dx is a symbol representing the difference in positions of
     * particle a and b.
     auto dx = create_dx(a, b);
     */
    
    
    /*
     * sum is a symbolic function that sums a sequence of 2d vectors
     AccumulateWithinDistance<std::plus<vdouble3>> sum;
     
     v[a] += dt * (
     //              // spring force between particles
     sum(b, if_else(id_[a] != id_[b],
     -prop->k_elast * (1.0 / norm(dx) - 1), 0.0) *
     dx));
     p[a] += dt * v[a];
     */
    
    //vdouble3 force;
    //vdouble3 torque;
    vdouble3 posj;
    vdouble3 orsj;
    vdouble3 posi;
    vdouble3 orsi;
    vdouble3 dxij;
    //vdouble3 ndx;
    vdouble3 meij;
    double nsqrij;
    int idi;
    int idj;
    double L=prop->L;
    NEIGHBOURS neis;
    //double dd_flat;
    //double dd_bend;
    double cc_flat;

    particles.init_neighbour_search(vdouble3(0,0,0),vdouble3(L,L,L),vbool3(false,false,false),prop->Rmax);
    particles.init_id_search();
    for (int i = 0; i < number; ++i) {
        idi=get<id>(particles[i]);
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        
        
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
            const typename particle_type::value_type& j = std::get<0>(tpl);
            idj=get<id>(j);
            
            if (idi<idj) {
                // j is const so we get it not const...
                typename particle_type::value_type jj=j;
                posj=get<position>(j);
                orsj=get<orientation>(j);
                dxij=posj-posi;
                //ndx=dxij/dxij.norm();
                //meij=orsi+orsj;

                nsqrij=dxij.squaredNorm();
               
                // Alignement force
                cc_flat=(dxij.dot(orsi+orsj))/(dxij.norm()*pow(nsqrij,3.0));
                get<force>(particles[i])+=cc_flat*orsi;
                get<force>(jj)-=cc_flat*orsj;
                
                
                // Bendind torque
                meij=prop->k_bend*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                get<torque>(particles[i])-=meij;
                get<torque>(jj)+=meij;

                // Lennard-Jones force
                //meij=dxij*((prop->k_rep)/pow(nsqrij,3)-(prop->k_att))/(pow(nsqrij,3.5));
                meij=dxij*((prop->k_rep)/pow(nsqrij,3.0)-(prop->k_att))/(pow(nsqrij,3.5));
                get<force>(particles[i])+=meij;
                get<force>(jj)-=meij;
                /*float iPerpend=-r_0*r_0*k_bend;
                 f_attract=fAtt[POS];
                 nd1.repForce.x+=fx;nd1.repForce.y+=fy;nd1.repForce.z+=fz;
                 nd2.repForce.x-=fx;nd2.repForce.y-=fy;nd2.repForce.z-=fz;
                 NF=(f_attract+f_rep);
                 fx=NF*dir.x;fy=NF*dir.y;fz=NF*dir.z;
                 */
                
                
                
                
                }

        }
    }

    
    
    
    
}

/* JAVA CODE FROM LIMESEG @KIARU
 
 
 
 float r6=(float)(java.lang.Math.pow((double)r_0,6d));
 float dx, dy, dz, norme2;
 dx=nd2.pos.x-nd1.pos.x;
 dy=nd2.pos.y-nd1.pos.y;
 dz=nd2.pos.z-nd1.pos.z;
 norme2=(dx*dx+dy*dy+dz*dz);
 if (norme2<limitInteractAttract) {
 float NF, fx, fy, fz;
 float f_attract, f_rep;
 int POS;
 float norme;
 boolean sameCell=nd1.ct.idInt==nd2.ct.idInt;
 dir.x=dx;dir.y=dy;dir.z=dz;
 norme=(float)(java.lang.Math.sqrt(norme2));
 dir.x/=norme;dir.y/=norme;dir.z/=norme;
 POS=(int)(norme/r_0*(float)(NStepPerR0));
 f_rep=fRep[POS];
 fx=f_rep*dir.x;fy=f_rep*dir.y;fz=f_rep*dir.z;
 if (sameCell) {
 f_attract=fAtt[POS];
 nd1.repForce.x+=fx;nd1.repForce.y+=fy;nd1.repForce.z+=fz;
 nd2.repForce.x-=fx;nd2.repForce.y-=fy;nd2.repForce.z-=fz;
 NF=(f_attract+f_rep);
 fx=NF*dir.x;fy=NF*dir.y;fz=NF*dir.z;
 
 nd1.N_Neighbor++;
 nd2.N_Neighbor++;
 // They are neighbor : we should align them compute their force
 // Computation of a force that tend to make a U shape between  Norm1/dir/Norm2 (moment / k_align)
 // So a force along Normal to flatten it proportional to the prodscal
 
 if (norme2<limitInteractBending) {
 
 SumN.x=nd1.Norm.x+nd2.Norm.x;SumN.y=nd1.Norm.y+nd2.Norm.y;SumN.z=nd1.Norm.z+nd2.Norm.z;
 float iFlatten=k_align*Vecteur3D.prodScal(SumN, dir)*(r6/((float)(java.lang.Math.pow((double)norme2,3d))));
 
 nd1.force.x+=nd1.Norm.x*iFlatten;
 nd1.force.y+=nd1.Norm.y*iFlatten;
 nd1.force.z+=nd1.Norm.z*iFlatten;
 
 nd2.force.x-=nd2.Norm.x*iFlatten;
 nd2.force.y-=nd2.Norm.y*iFlatten;
 nd2.force.z-=nd2.Norm.z*iFlatten;
 // We compute the perpendicular with the vector dir
 float iPerpend=-r_0*r_0*k_bend;
 Vecteur3D couple=Vecteur3D.prodVect(nd1.Norm,nd2.Norm);
 couple.scale(r6/((float)(java.lang.Math.pow((double)norme2,3d))));
 //Vecteur3D.prodScal(nd1.Norm,dir);
 nd1.moment.x+=iPerpend*couple.x;
 nd1.moment.y+=iPerpend*couple.y;
 nd1.moment.z+=iPerpend*couple.z;
 //.float iPerpend2=-r_0*r_0*k_bend*Vecteur3D.prodScal(nd2.Norm,dir);
 nd2.moment.x-=iPerpend*couple.x;
 nd2.moment.y-=iPerpend*couple.y;
 nd2.moment.z-=iPerpend*couple.z;
 }
 
 float iPerpend1=-r_0*r_0*k_bend*Vecteur3D.prodScal(nd1.Norm,dir);
 nd1.moment.x+=iPerpend1*dir.x;
 nd1.moment.y+=iPerpend1*dir.y;
 nd1.moment.z+=iPerpend1*dir.z;
 float iPerpend2=-r_0*r_0*k_bend*Vecteur3D.prodScal(nd2.Norm,dir);
 nd2.moment.x+=iPerpend2*dir.x;
 nd2.moment.y+=iPerpend2*dir.y;
 nd2.moment.z+=iPerpend2*dir.z;
 */





void Part_set::NextStep(const Meshless_props* simul_prop){
    float dt=simul_prop->dt;
    
    
    
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