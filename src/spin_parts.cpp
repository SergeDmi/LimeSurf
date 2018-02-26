#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include <iostream>
#include "assert_macro.h"
//using namespace tinyply;
#include "tinyply.h"
using namespace tinyply;


const double PI = boost::math::constants::pi<double>();

Part_set::Part_set(Part_set_props * p) {
    // Dummy creator
    number=0;
    prop=p;
    double L=prop->L;
    particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
}

void Part_set::create() {
    
    if (prop->init_shape==0){
        number=PutOnSphere();
    }
    
    if (prop->init_shape==1){
        number=PutOnSheet();
    }
    
     std::cout << "# created " << number << "particles, with expected R0 " << prop->R0 << std::endl;
}



void Part_set::create(std::string fname) {
    
    number=load_from_text(fname);
    
    std::cout << "# created " << number << "particles, with expected R0 " << prop->R0 << std::endl;
}


int Part_set::load_from_text(std::string fname){
    //std::string fname;
    //fname = source;
    std::ifstream ss(fname, std::ios::binary);
    if (ss.fail())
    {
        throw std::runtime_error("failed to open " + fname);
    }
    
    PlyFile file;
    file.parse_header(ss);
    std::shared_ptr<PlyData> vertices, normals, colors, faces, texcoords;
    
    try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
    
    try { normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
    
    file.read(ss);
    int nv=vertices->count;
    
    if (normals->count != nv) {
        std::cerr << " Error : number of vertices and normals should be the same" << std::endl;
    }
    else
    {
            std::cout << "# \tRead " << vertices->count << " total vertices "<< std::endl;
    }
    
    vdouble3 zero(0,0,0);
    struct float3 { float x, y, z; };

    std::vector<float3> verts(nv);
    std::vector<float3> norms(normals->count);
    const size_t numVerticesBytes = vertices->buffer.size_bytes();
    const size_t numNormals_Bytes =  normals->buffer.size_bytes();
    std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
    std::memcpy(norms.data(),  normals->buffer.get(), numNormals_Bytes);
    
    //std::cout << "compte : " << nv << "numVerticesBytes "  << numVerticesBytes << std::endl;
    
    //std::cout << "try : " << vertices->t << std::endl;
    //if (vertices->t == tinyply::Type::FLOAT32) { std::cout << "FLOAT" << std::endl; }
    //if (vertices->t == tinyply::Type::FLOAT64) { std::cout << "DOUBLE" << std::endl; }
    
    
    //std::cout << verts[0].y << " " << std::endl;
    //std::memcpy(norms.data(),  normals->buffer.get(), numNormals_Bytes);
    
    
    for (int i=0; i<nv; ++i) {
        typename particle_type::value_type p;
        vdouble3 pos(verts[i].x,verts[i].y,verts[i].z);
        vdouble3 dir(norms[i].x,norms[i].y,norms[i].z);
        
        //std::cout << dir << " " << std::endl;
        
        get<position>(p) = pos;
        get<orientation>(p) = dir;
        get<force>(p) = zero;
        get<torque>(p) = zero;
        
        particles.push_back(p);

    }
    

    
    return nv;
}



int Part_set::num() {return number;}


// Put many particles on a sphere at random
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

    for ( i = 0; i < N; ++i) {
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

// Put many particles on a sheet with hexagonal lattice
int Part_set::PutOnSheet(){
    int N=prop->init_number;
    double R=prop->init_radius;
    
    double L=prop->L;
    vdouble3 pos1;
    vdouble3 pos2;
    vdouble3 up(0,0,1);
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
            pos1=vdouble3(i*h*2,j*R0,2*R0*sin((2*PI*i*h*2)/R));
            pos2=vdouble3(i*h*2+h,j*R0+0.5*R0,2*R0*sin((2*PI*(i*h*2+h))/R));
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
    
   
    //number=count;
    return count;
}

// Makes sure everything is in place
void Part_set::GetStarted(){
    double L=prop->L;
    particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
    //particles.update_positions();
    if (prop->elastic) {
        particles.init_id_search();
        GetNeighbours();
    }
    
    //particles.init_id_search();
}

// Needed when elastic will be implemented
void Part_set::GetNeighbours() {
    int count=0;
    int ide;
    int n;
    for (int i = 0; i < number; ++i) {
        n=0;
        std::vector<int> neis;
        int idi=get<id>(particles[i]);
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
            const typename particle_type::value_type& j = std::get<0>(tpl);
            ide=get<id>(j);
            if (ide!=idi) {
                    neis.push_back(ide);
                    n++;
            }
        }

        get<nn>(particles[i])=n;
        get<neighbours>(particles[i])=neis;
        if (n>count) {
            count=n;
        }
            
    }
    //std::cout << "# last vec" << std::endl;
    //for(int k : neis) {
    //    std::cout << "#              " << k << std::endl;
    //}
    
    max_neighbours=count;
    std::cout << "#Maximum neighbours number : " << max_neighbours << std::endl;
 
}

// One simulation step
void Part_set::NextStep(const Meshless_props* simul_prop){
    ComputeForces();
    IntegrateForces(simul_prop);
    if (!(prop->elastic)) {
        particles.update_positions();
    }
    ClearForces();
}

void Part_set::ComputeForces(){
    if (prop->elastic) {
        ComputeForcesElastic();
    }
    else
    {
        ComputeForcesViscous();
    }
}

// Erazing forces
void Part_set::ClearForces() {
    vdouble3 zero(0,0,0);
    for (int i = 0; i < number; ++i) {
        get<force>(particles[i])=zero;
        get<torque>(particles[i])=zero;
    }
    
}

// Adding forces to particles
void Part_set::IntegrateForces(const Meshless_props* simul_prop){
    float dt=simul_prop->dt;
    vdouble3 forcei;
    for (int i = 0; i < number; ++i) {
        // Upper threshold for forces
        //      Ugly but needed if no minimum distance
        forcei=get<force>(particles[i]);
        if (forcei.squaredNorm()>prop->Fmax2) {
            forcei=forcei*((prop->Fmax)/forcei.norm());
            //std::cout << "# warning : over the top ; now : " << forcei << std::endl;
        }
        get<position>(particles[i])+=((forcei)*(dt/prop->visco));
        //std::cout << "# checking forcei : " << forcei << std::endl;
        get<orientation>(particles[i])+=cross(get<orientation>(particles[i]),get<torque>(particles[i]))*(dt/prop->Rvisc);
    }
    
}


void Part_set::ComputeForcesElastic(){
    vdouble3 posi;
    vdouble3 orsi;
    int idi,idj;
    double k_elast=prop->k_elast;
    double norm;
    vdouble3 dir;
    vdouble3 tri;
    double proj;
    int count;
    for (int i = 0; i < number; ++i) {
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        idi=get<id>(particles[i]);
        int neibs=get<nn>(particles[i]);
        std::vector<int> neis=get<neighbours>(particles[i]);
        std::random_shuffle ( neis.begin(), neis.end() );
        count=0;
        int R2mean=0;
        vdouble3 oldvec(0,0,0);
        //std::cout << "# neibs=" << neibs << std::endl;
        //for (int j=0; j<neibs; ++j) {
        //for (std::vector<int>::iterator it=neis.begin(); it!=neis.end(); ++it) {
        //    idj=*it;
        for(int idj : neis) {
        //for (int bj=0; bj<neibs; ++bj) {
        //    idj=neis[bj];
            auto j = particles.get_query().find(idj);
            //auto id_2_ref = *particles.get_query().find(2);
            //auto id_2_position_reference = *get<position>(id_2);
            vdouble3 posj=*get<position>(j);
            //vdouble3 orsj=*get<orientation>(idj);
            //vdouble3 sumo=orsi+orsj;
            vdouble3 dxij=posj-posi;
            norm=dxij.norm();
            dir=dxij*prop->R0/norm;
            get<force>(particles[i])+=dxij-dir;
            //get<force>(particles[i])+=dxij-dir;
            
            R2mean+=norm*norm;
            // This is highly suspicious !
            // We compute the normals from triangles
            // To update the normal of the point
            // Is it any good ?
            if (count>0) {
                tri=cross(dir,oldvec);
                proj=orsi.dot(tri);
                tri*=proj/tri.norm();
                 // Problematic :
                get<torque>(particles[i])-=tri;
            }
            oldvec=dir;
            count++;
        }
        // Problematic :
        //get<force>(particles[i])+=orsi*(prop->pressure*R2mean)/count;
        
    }
    //std::cout << "# computed " << count << " neighbours" <<std::endl;
}



// Computing forces with viscous setting (i.e. changeable nearest neighbours)
void Part_set::ComputeForcesViscous(){
    vdouble3 posi;
    vdouble3 orsi;
    int idi,idj;
    double L=prop->L;
    //NEIGHBOURS neis;
    double cc_flat;
    double p_att=(1.0+prop->p_att)/2.0;
    double p_align=(1.0+prop->p_align)/2.0;
    double p_rep=(prop->p_rep+prop->p_att)/2.0;
    int neibs;
    for (int i = 0; i < number; ++i) {
        neibs=0;
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        idi=get<id>(particles[i]);
        for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
             typename particle_type::value_type j = std::get<0>(tpl);
            
            idj=get<id>(j);
            if (idi!=idj) {
                neibs++;
                vdouble3 posj=get<position>(j);
                vdouble3 orsj=get<orientation>(j);
                vdouble3 sumo=orsi+orsj;
                vdouble3 dxij=posj-posi;
                double nsqrij=std::max(dxij.squaredNorm(),prop->minR);
                // Bending torque
                get<torque>(particles[i])-=(prop->k_bend)*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                // Forces : Lennard Jones & alignement
                get<force>(particles[i])+=(-dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att))
                                       +sumo*(0.5*prop->k_align*(dxij.dot(sumo)/pow(nsqrij,p_align))));
            //assert_true(dxij.dot(sumo)>0.0);
            //get<force>(particles[i])+=sumo*(0.5*prop->k_align*(dxij.dot(sumo)/pow(nsqrij,p_align)));
            }
        }
        
        get<nn>(particles[i])=neibs;
        //std::cout << "# neibs=" << neibs << std::endl;
    }
}


// Temporary dirty exporting to a file
void Part_set::Export(int t){
#ifdef HAVE_VTK
    vtkWriteGrid("particles",t,particles.get_grid(true));
#endif
    ofstream exportfile;
    exportfile.open ("particles_out.txt");
    //myfile << "Writing this to a file.\n";
    exportfile << "# X  Y   Z   Fx  Fy  Fz Nneigh ID \n";
    for (int i = 0; i < number; ++i) {
        exportfile <<" "<< get<position>(particles[i])[0] <<" "<< get<position>(particles[i])[1] << " "<< get<position>(particles[i])[2] << " " << get<orientation>(particles[i])[0] << " " << get<orientation>(particles[i])[1] << " " << get<orientation>(particles[i])[2] <<" "<< get<nn>(particles[i])<<" " << get<id>(particles[i])<< " \n";
        
    }
    exportfile.close();
}

