#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include <iostream>
#include <fstream>
#include "assert_macro.h"
#include <limits>
#include "tinyply.h"
#include "spin_parts_props.h"
using namespace tinyply;
using namespace std;
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep


const double PI = boost::math::constants::pi<double>();

// Dummy constructor
Part_set::Part_set() {
    number=0;
    diverging=false;
   
}   

Part_set::Part_set(Part_set_props * p) {
    number=0;
    prop=p;
    double L=prop->L;
    diverging=false;
   
    //double x_conf=prop->x_conf;
    //double y_conf=prop->y_conf;
    //double z_conf=prop->z_conf;
    // Makin da functions yo
    /*
    if (prop->x_conf>0) {
        std::cout << "# There will be confinement along x" << std::endl;
        add_x_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[0]-=pp->x_conf*((posi[0]>pp->x_max)*(posi[0]-pp->x_max) + (posi[0]<-pp->x_max)*(posi[0]+pp->x_max)); };
    } else {
        add_x_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    
    if (prop->y_conf>0) {
        std::cout << "# There will be confinement along y" << std::endl;
        add_y_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[1]-=pp->y_conf*((posi[1]>pp->y_max)*(posi[1]-pp->y_max) + (posi[1]<-pp->y_max)*(posi[1]+pp->y_max)); };
    } else {
        add_y_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    
    if (prop->z_conf>0) {
        std::cout << "# There will be confinement along z" << std::endl;
        add_z_conf = [pp=prop] (vdouble3& force, vdouble3 posi ) { force[2]-=pp->z_conf*((posi[2]>pp->z_max)*(posi[2]-pp->z_max) + (posi[2]<-pp->z_max)*(posi[2]+pp->z_max)); };
    } else {
        add_z_conf = [] (vdouble3& force, vdouble3 posi ) {};
    }
    */
}

// If particles need to be created from properties
void Part_set::create() {
    std::cout << "# starting particle creation ; load " << prop->fname_in << " if " << prop->load_from_file << std::endl;
    if (prop->load_from_file>0) {
        std::string test=prop->fname_in;
        std::cout << "# ... Trying to load " << test << std::endl;
        number=load_from_text();
        
    }
    else {
    if (prop->init_shape==0){
        number=PutOnSphere();
    }
    
    if (prop->init_shape==1){
        number=PutOnSheet();
    }
    }
     std::cout << "# created " << number << "particles, with expected R0 " << prop->R0 << std::endl;
     
}



// Loading from ply file
// Copied from tinyply's example
// @TODO : make sure that memcpy is cosher
int Part_set::load_from_text(){
    //std::cout << "# ... acquiring file " << prop->fname_in << std::endl;
    std::ifstream ss(prop->fname_in, std::ios::binary);
    if (ss.fail())
    {
        throw std::runtime_error("failed to open " + prop->fname_in);
    }
    PlyFile file;
    file.parse_header(ss);
    std::shared_ptr<PlyData> vertices, normals, colors, faces, tetra,texcoords;
    bool is_tetra=1;
    
    try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
    
    try { normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
    
    try { faces = file.request_properties_from_element("face", { "vertex_index" }); }
    catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }
    
    try { tetra = file.request_properties_from_element("tetrahedra", { "vertex_index" }); }
    catch (const std::exception & e) { is_tetra=0; }
    
    file.read(ss);
    int nv=vertices->count;
    int nf=faces->count;
    int nt;
    if (is_tetra) {
        nt=tetra->count;
    }
    else {nt=0;}
    
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
    
    const size_t numFacesBytes = faces->buffer.size_bytes();
    face_list tris(nf);
    std::memcpy(tris.data(), faces->buffer.get(), numFacesBytes);
    triangles=tris;
    
    tetr_list quatrs(nt);
    if (is_tetra) {
        const size_t numTetrasBytes = tetra->buffer.size_bytes();
        std::memcpy(quatrs.data(), tetra->buffer.get(), numTetrasBytes);
        
    }
    tetrahedra=quatrs;
    
    const size_t numVerticesBytes = vertices->buffer.size_bytes();
    const size_t numNormals_Bytes =  normals->buffer.size_bytes();
    std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
    std::memcpy(norms.data(),  normals->buffer.get(), numNormals_Bytes);
    
    
    for (int i=0; i<nv; ++i) {
        typename particle_type::value_type p;
        vdouble3 pos(verts[i].x,verts[i].y,verts[i].z);
        vdouble3 dir(norms[i].x,norms[i].y,norms[i].z);
        
        get<position>(p) = pos;
        get<orientation>(p) = dir;
        get<force>(p) = zero;
        get<torque>(p) = zero;
        get<state>(p) = -1.0;
        get<nn>(p) = 0;
        particles.push_back(p);
        
    }
    
    n_faces=nf;
    n_tetras=nt;
    
    return nv;
}

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
            
            // Using aboria's power
            
            for (auto tpl=euclidean_search(particles.get_query(),get<position>(p),mindist); tpl!=false ; ++tpl) {
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
    return count;
}
    

// Makes sure everything is in place
void Part_set::GetStarted(){
    CheckBoxSize();
    particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
    
}

// Check if box size is OK
// @TODO : repair !
void Part_set::CheckBoxSize() {
    vdouble3 bottomleft(INFINITY,INFINITY,INFINITY);
    vdouble3 topeuright(-INFINITY,-INFINITY,-INFINITY);
    for (int i = 0; i < number; ++i) {
        vdouble3 posi=get<position>(particles[i]);
        for (int ix=0; ix<3; ++ix) {
            if (posi[ix]<bottomleft[ix]) {
                bottomleft[ix]=posi[ix];
            }
            else if (posi[ix]>topeuright[ix]) {
                topeuright[ix]=posi[ix];
              
            }
        }
    }
    
    for (int ix=0; ix<3; ++ix) {
        double dix=topeuright[ix]-bottomleft[ix];
        bottomleft[ix]=bottomleft[ix]-0.15*dix;
        topeuright[ix]=topeuright[ix]+0.15*dix;
        if (bottomleft[ix]<prop->corner_0[ix]) { prop->corner_0[ix]=bottomleft[ix];}
        if (topeuright[ix]>prop->corner_1[ix]) { prop->corner_1[ix]=topeuright[ix];}
    }
        std::cout << "# bounding box from " << prop->corner_0 << " to " << prop->corner_1 << std::endl;
}

// Finds the furthest points in x,x, y,y, z,z
void Part_set::FindBounds() {
    vdouble3 bottomleft(INFINITY,INFINITY,INFINITY);
    vdouble3 topeuright(-INFINITY,-INFINITY,-INFINITY);
    // Looping on all points to find the min & max of x y z
    for (int i = 0; i < number; ++i) {
        vdouble3 posi=get<position>(particles[i]);
        for (int ix=0; ix<3; ++ix) {
            if (posi[ix]<bottomleft[ix]) {
                bottomleft[ix]=posi[ix];
            }
            else if (posi[ix]>topeuright[ix]) {
                topeuright[ix]=posi[ix];   
            }
        }
    }
    // Now we can update bounds
    for (int ix=0; ix<3; ++ix) {
        bounds[ix]=bottomleft[ix];
        bounds[3+ix]=topeuright[ix];
    }
}

// Get the closest neighbours
void Part_set::GetNeighbours() {
    int count=0;
    int tnuoc=666;
    int ide;
    int n;
    vdouble3 dir(0,0,0);
    vdouble3 pos1(0,0,0);
    
    for (int i = 0; i < number; ++i) {
        n=0;
        neigh_pairs pairs;
        pos1=get<position>(particles[i]);
        int idi=get<id>(particles[i]);
        //for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
        for (auto tpl = euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax); tpl!=false ; ++tpl ) {
            //const typename particle_type::value_type& j = std::get<0>(tpl);
            //const auto j = std::get<0>(*tpl);
            ide=get<id>(*tpl);
            if (ide!=idi) {
                    dir=get<position>(*tpl)-pos1;
                    pair_n ppp(ide,sqrt(dir.squaredNorm()));
                    pairs.push_back(ppp);
                    n++;
            }
        }

        while (n>prop->max_neighbours) {
            n=pop_furthest_neighbour(&pairs,n);
        }
        get<nn>(particles[i])=n;
        get<neighbours>(particles[i])=pairs;
        
        if (n>count) {
            count=n;
        }
        if (n<tnuoc) {
            tnuoc=n;
        }
            
    }
    
    max_neighbours=count;
  
}

// we should check we got a decent part set after loading from text
void Part_set::CheckPartSet() {
    int count=0;
    int tnuoc=666;
    vdouble3 orsi;
    vdouble3 zero(0,0,0);
    double normi;
    for (int i = 0; i < number; ++i) {
        orsi=get<orientation>(particles[i]);
        normi=orsi.squaredNorm();
        if (std::isnan(normi)) {
            get<orientation>(particles[i])=zero;
        }
        neigh_pairs pairs=get<neighbours>(particles[i]);
        int n=static_cast<int>(pairs.size());
        
        if (n>count) {
            count=n;
        }
        if (n<tnuoc) {
            tnuoc=n;
        }
    }
    std::cout << "# Maximum neighbours number : " << count << std::endl;
    std::cout << "# Minimum neighbours number : " << tnuoc << std::endl;
}

// removes furthest neighbour from particle
int Part_set::pop_furthest_neighbour(neigh_pairs * pairs, int n) {
    double dist=0;
    int ex=-1;
    for (int i=0;i<n;++i)
    {
        pair_n ppp=pairs->at(i);
        if (ppp.second>dist) {
            dist=ppp.second;
            ex=i;
        }
    }
    pairs->erase(pairs->begin()+ex);
    return n-1;
}


// One simulation step
void Part_set::NextStep(const Simul_props & simul_prop){
    ComputeForces(simul_prop);
    AddConfinementForces(simul_prop);
    IntegrateForces(simul_prop);
    if (!(prop->elastic)) {
        particles.update_positions();
    }
    ClearForces();
}


// Setting forces & torques to 0
void Part_set::ClearForces() {
    vdouble3 zero(0,0,0);
    for (int i = 0; i < number; ++i) {
        get<force>(particles[i])=zero;
        get<torque>(particles[i])=zero;
    }
    
}

// Adding forces to particles
void Part_set::IntegrateForces(const Simul_props & simul_prop){
    double dt_trans=simul_prop.dt/prop->visco;
    double dt_rot=simul_prop.dt/prop->Rvisc;
    vdouble3 forcei;
    for (int i = 0; i < number; ++i) {
        
        forcei=get<force>(particles[i]);
        
        // We didn't implement maximum force.
        //if (forcei.squaredNorm()>prop->Fmax2) {
        //    forcei=forcei*((prop->Fmax)/forcei.norm());
            //std::cout << "# warning : over the top ; now : " << forcei << std::endl;
        //}
        
        // Applying force & torque
        get<position>(particles[i])+=((forcei)*(dt_trans));
        get<orientation>(particles[i])+=cross(get<orientation>(particles[i]),get<torque>(particles[i]))*(dt_rot);
    }
}

// Renormalizing the normals
void Part_set::RenormNorms(){
    vdouble3 normi;
    for (int i = 0; i < number; ++i) {
        normi=get<orientation>(particles[i]);
        if ((get<state>(particles[i])>0)) {
            get<orientation>(particles[i])/=sqrt(normi.squaredNorm());            
        }
    }
}

// Computing forces with viscous setting (i.e. changeable nearest neighbours)
//@TODO : verify this part
void Part_set::ComputeForces(const Simul_props & simul_prop){
    
    vdouble3 posi;
    vdouble3 orsi;
    int idi,idj;
    double p_att=(1.0+prop->p_att)/2.0;
    double p_align=(1.0+prop->p_align)/2.0;
    double p_rep=(prop->p_rep+prop->p_att)/2.0;
    int neibs;
    for (int i = 0; i < number; ++i) {
        neibs=0;
        posi=get<position>(particles[i]);
        orsi=get<orientation>(particles[i]);
        idi=get<id>(particles[i]);
        //for (auto tpl: euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax)) {
        for (auto tpl = euclidean_search(particles.get_query(),get<position>(particles[i]),prop->Rmax); tpl!=false ; ++tpl ) {
             //typename particle_type::value_type j = std::get<0>(*tpl);
             //const typename particle_type::value_type& j = std::get<0>(*tpl);
            
            idj=get<id>(*tpl);
            if (idi!=idj) {
                neibs++;
                vdouble3 posj=get<position>(*tpl);
                vdouble3 orsj=get<orientation>(*tpl);
                vdouble3 sumo=orsi+orsj;
                vdouble3 dxij=posj-posi;
                double nsqrij=std::max(dxij.squaredNorm(),prop->minR);
                // Bending torque
                get<torque>(particles[i])-=(prop->k_bend)*cross(orsi,orsj)/(pow(nsqrij,3.0));;
                // Forces : Lennard Jones & alignement
                get<force>(particles[i])+=(-dxij*((prop->k_rep)/pow(nsqrij,p_rep)-(prop->k_att))/(pow(nsqrij,p_att))
                                       +sumo*(0.5*prop->k_align*(dxij.dot(sumo)/pow(nsqrij,p_align))));
            }
        }
        
        get<nn>(particles[i])=neibs;
    }
     //* */
}

// Confinement forces
void Part_set::AddConfinementForces(const Simul_props & simul_prop){
    vdouble3 posi;
    //vdouble3 orsi;
    vdouble3 forcei;
    vdouble3 zero(0,0,0);
    //int idi,idj;
    //double p_att=(1.0+prop->p_att)/2.0;
    //double p_align=(1.0+prop->p_align)/2.0;
    //double p_rep=(prop->p_rep+prop->p_att)/2.0;
    //int neibs;
    
    //if (prop->x_conf+prop->y_conf+prop->z_conf>0) {
    if ( simul_prop.is_confinement ) {        
        // There is actually confinement
        //for (int i = 0; i < number; ++i) {
        for (auto part : particles) {
            forcei=zero;
            posi=get<position>(part);
            simul_prop.add_confinement_force(forcei,posi);
            //add_x_conf(forcei,posi);
            //add_y_conf(forcei,posi);
            
            //std::cout << " force bla "  << forcei[0] << " , " << forcei[1] << " , " << forcei[2] << std::endl;
            //add_z_conf(forcei,posi);
            //idi=get<id>(particles[i]);
            
            get<force>(part)+=forcei;
            
            //std::cout << " force i "  << forcei[0] << " , " << forcei[1] << " , " << forcei[2] << std::endl;
        }
    }
    
}



// Temporary dirty exporting to a file
void Part_set::Export(int t){
#ifdef HAVE_VTK
    //vtkWriteGrid("particles",t,particles.get_grid(true));
#endif
    ofstream exportfile;
    std::string fname_out="particles_out_"+std::to_string(t)+".txt";
    exportfile.open(fname_out);
    //myfile << "Writing this to a file.\n";
    exportfile << "# X  Y Z   Nx  Ny  Nz Fx Fy Fz Nneigh ID \n";
    for (int i = 0; i < number; ++i) {
        exportfile <<" "<< get<position>(particles[i])[0] <<" "<< get<position>(particles[i])[1] << " "<< get<position>(particles[i])[2];
        exportfile<< " " << get<orientation>(particles[i])[0] << " " << get<orientation>(particles[i])[1] << " " << get<orientation>(particles[i])[2];
        exportfile<< " " << get<force>(particles[i])[0] << " " << get<force>(particles[i])[1] << " " << get<force>(particles[i])[2];
        exportfile <<" "<< get<nn>(particles[i])<<" " << get<id>(particles[i])<< " \n";
        
    }
    exportfile.close();
}

// saving to .ply
// @TODO : problem in the tags of exported ply files
void Part_set::Export_bly(int n_frame,const Simul_props & simul_prop, double t){
    std::string fname=prop->fname_out+simul_prop.name;
    //std::ifstream ss(fname, std::ios::binary);
    std::string filename;
    std::filebuf fb;
    std::string numero(std::to_string(n_frame));
    std::string uno(std::to_string(simul_prop.n_frames));
    while (numero.length()<uno.length()) {
        numero="0"+numero;
    }
    filename=fname+"_"+numero+".ply";
    fb.open(filename, std::ios::out | std::ios::binary);
    std::ostream outputStream(&fb);
    PlyFile exampleOutFile;
    vdouble3 orsi;
    vdouble3 posi;
   
    // First we fill in the vertices and normals
    struct float3 { float x, y, z; };
    std::vector<float3> verts(number);
    std::vector<float3> norms(number);
    //std::cout << "number of vertices : " << number << std::endl;
    for (int idi = 0; idi < number; ++idi) {
        // We search by ID to make sure no particle reshuffling occured
        auto i = particles.get_query().find(idi);
        posi=*get<position>(i);
        orsi=*get<orientation>(i);
        //verts[idi]={(float)posi[0],(float)posi[1],(float)posi[2]};
        verts[idi].x=(float)posi[0];verts[idi].y=(float)posi[1];verts[idi].z=(float)posi[2];
        norms[idi].x=(float)orsi[0];norms[idi].y=(float)orsi[1];norms[idi].z=(float)orsi[2];
        //norms[idi]={(float)orsi[0],(float)orsi[1],(float)orsi[2]};
        //std::cout << " "  << idi << std::flush;
    }
    
    exampleOutFile.add_properties_to_element("vertex", { "x", "y", "z" }, 
                            Type::FLOAT32, verts.size(), reinterpret_cast<uint8_t*>(verts.data()), Type::INVALID, 0);
    exampleOutFile.add_properties_to_element("vertex", { "nx", "ny", "nz" }, 
                            Type::FLOAT32, verts.size(), reinterpret_cast<uint8_t*>(norms.data()), Type::INVALID, 0);
    exampleOutFile.add_properties_to_element("face", { "vertex_index" },
                            Type::UINT32,  triangles.size(), reinterpret_cast<uint8_t*>(triangles.data()), Type::UINT8, 3);
                            
    exampleOutFile.get_comments().push_back("generated by tinyply");
    exampleOutFile.get_comments().push_back("time "+std::to_string(t));
    exampleOutFile.write(outputStream, false);
        
    fb.close();

}

