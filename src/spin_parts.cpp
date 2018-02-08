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



void Part_set::NextStep(Meshless_props* simul_prop){
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