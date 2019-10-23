#include <random>
#include "Aboria.h"
using namespace Aboria;
#include <boost/math/constants/constants.hpp>
#include <math.h>
#include "spin_parts.h"
#include <iostream>
#include <fstream>
//#include "assert_macro.h"
#include <limits>
#include "tinyply.h"
#include "spin_parts_props.h"
using namespace tinyply;
using namespace std;
#include "yaml-cpp/yaml.h" 

/*
 A class storing a set of particles ; basis for all wall representations
*/

// Not even sure we need this here
const double PI = boost::math::constants::pi<double>();

// Dummy constructor
Part_set::Part_set() {
    number=0;
    diverging=false;
}   

// Constructor from properties
Part_set::Part_set(Part_set_props * p) {
    number=0;
    prop=p;
    diverging=false;
}

// If particles need to be created from properties
void Part_set::create() {
    // std::cout << "# starting particle creation ; load " << prop->fname_in << std::endl;
    // std::string test=prop->fname_in;
    // std::cout << "# ... Trying to load " << test << std::endl;
    number=load_from_text();
    std::cout << "# created " << number << "particles" << std::endl;
     
}


// Finds the furthest points in x,x, y,y, z,z
// @TODO : needs to be fixed
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

// Loading from ply file
// @TODO : make sure that memcpy is cosher
int Part_set::load_from_text(){
    std::cout << "# ... acquiring file " << prop->fname_in << std::endl;
    // Readint the desired file
    std::ifstream ss(prop->fname_in, std::ios::binary);
    if (ss.fail())
    {
        // Well that doesn't start well now, does it ?
        throw std::runtime_error("failed to open " + prop->fname_in);
    }
    
    // Copied from tinyply's example
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
        std::cout << "# \tRead " << tetra->count << " total tetrahedrons "<< std::endl;
    }
    else {nt=0;}
    
    if (normals->count != nv) {
        std::cerr << " Error : number of vertices and normals should be the same" << std::endl;
    }
    else
    {
            std::cout << "# \tRead " << vertices->count << " total vertices "<< std::endl;
    }

    // Let us not create zero too often
    vdouble3 zero(0,0,0);
    
    // a struct containint 3 floats
    struct float3 { float x, y, z; };
    // and vectors thereof
    std::vector<float3> verts(nv);
    std::vector<float3> norms(normals->count);
    // we will fill the data of these vectors from the buffer
    // face list is a vector of structures of 3 ints
    face_list tris(nf);
    // memory required
    const size_t numFacesBytes = faces->buffer.size_bytes();
    // Dirty copy, memcpy
    std::memcpy(tris.data(), faces->buffer.get(), numFacesBytes);
    // But we got our triangles
    triangles=tris;
    
    // Reading tetrahedrons if we do have any
    tetr_list quatrs(nt);
    if (is_tetra) {
        const size_t numTetrasBytes = tetra->buffer.size_bytes();
        std::memcpy(quatrs.data(), tetra->buffer.get(), numTetrasBytes);
    }
    tetrahedra=quatrs;
    
    // Copying memory buffers for vertices and normals
    const size_t numVerticesBytes = vertices->buffer.size_bytes();
    const size_t numNormals_Bytes =  normals->buffer.size_bytes();
    std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
    std::memcpy(norms.data(),  normals->buffer.get(), numNormals_Bytes);
    
    // Now we actually make the particles !
    for (int i=0; i<nv; ++i) {
        typename particle_type::value_type p;
        vdouble3 pos(verts[i].x,verts[i].y,verts[i].z);
        vdouble3 dir(norms[i].x,norms[i].y,norms[i].z);
     
        // p is a particle, let's fill it
        get<position>(p) = pos;
        get<orientation>(p) = dir;
        get<force>(p) = zero;
        get<torque>(p) = zero;
        get<state>(p) = -1.0;
        get<nn>(p) = 0;
        // now we add it to our vector of particles
        particles.push_back(p);
        
        // Let's keep this here, you know, just in cae
        /* 
        if (i==0) {
            std::cout << "Position     :  " << pos[0] << " , " << pos[1] << " , " << pos[2] << std::endl;
            std::cout << "Orientation  :  " << dir[0] << " , " << dir[1] << " , " << dir[2] << std::endl;
        }
         */
    }
   
    // Book keeping
    n_faces=nf;
    n_tetras=nt;

    // Done !
    return nv;
}
    

// Makes sure everything is in place
// @TODO : make it actually do something
void Part_set::GetStarted(){
    // We shoudl do an init_neighbour search if we plan to neighbour search.
    // Removed because there are problems with checkboxsize
    //CheckBoxSize();
    //particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    //std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
}



// we should check we got a decent part set after loading from text
void Part_set::CheckPartSet() {
    int count=0;
    int tnuoc=INT_MAX;
    vdouble3 orsi;
    vdouble3 zero(0,0,0);
    double normi;
    int n;
    for (int i = 0; i < number; ++i) {
        // Trying to get a non-isnan orientation !
        orsi=get<orientation>(particles[i]);
        normi=orsi.squaredNorm();
        if (std::isnan(normi)) {
            get<orientation>(particles[i])=zero;
        }
        // Counting min and max number of neighbours
        n=get<nn>(particles[i]);
        
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


// One simulation step
// @TODO : make {} ?
void Part_set::NextStep(const Simul_props & simul_prop){
    // This should be virtualized !  
    ComputeForces(simul_prop);
    AddConfinementForces(simul_prop);
    IntegrateForces(simul_prop);
    ClearForces();
}


// Setting forces & torques to 0
void Part_set::ClearForces() {
    vdouble3 zero(0,0,0);
    for (auto part : particles) {
        get<force> (part)=zero;
        get<torque>(part)=zero;
    }   
}

// Adding forces to particles
void Part_set::IntegrateForces(const Simul_props & simul_prop){
    double dt_trans=simul_prop.dt/prop->visco;
    double dt_rot=simul_prop.dt/prop->Rvisc;

    for (auto part : particles) {        
        // Applying force & torque
        get<position>(part)   +=(get<force>(part)*(dt_trans));
        get<orientation>(part)+=cross(get<orientation>(part),get<torque>(part))*(dt_rot);
    }
}

// Renormalizing the normals
void Part_set::RenormNorms(){
    for (auto part : particles) {  
        if ((get<state>(part)>0)) {
            get<orientation>(part)/=sqrt(get<orientation>(part).squaredNorm());            
        }
    }
}

// Confinement forces
void Part_set::AddConfinementForces(const Simul_props & simul_prop){
    vdouble3 posi;
    vdouble3 forcei;
    vdouble3 zero(0,0,0);
    
    if ( simul_prop.is_confinement ) {        
        // There is actually confinement
        for (auto part : particles) {
            forcei=zero;
            posi=get<position>(part);
            simul_prop.add_confinement_force(forcei,posi);
            
            get<force>(part)+=forcei;
          
        }
    }
    
}



// Dirty exporting to a file
// Kept for safety
void Part_set::Export(int t, double time){
#ifdef HAVE_VTK
    vtkWriteGrid("particles",t,particles.get_grid(true));
#endif
    ofstream exportfile;
    std::string fname_out="particles_out_"+std::to_string(t)+".txt";
    exportfile.open(fname_out);
    exportfile << "# time : " << time << std::endl;
    exportfile << "# visco : " << prop->visco << std::endl;
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
    exampleOutFile.get_comments().push_back("visco "+std::to_string(prop->visco));
    exampleOutFile.get_comments().push_back("pressure "+std::to_string(simul_prop.pressure));
    exampleOutFile.write(outputStream, false);
        
    fb.close();

}

