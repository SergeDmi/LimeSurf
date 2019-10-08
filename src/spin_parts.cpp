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
    //double L=prop->L;
    diverging=false;
  
}

// If particles need to be created from properties
void Part_set::create() {
    std::cout << "# starting particle creation ; load " << prop->fname_in << std::endl;
        std::string test=prop->fname_in;
        std::cout << "# ... Trying to load " << test << std::endl;
        number=load_from_text();
    
     std::cout << "# created " << number << "particles" << std::endl;
     
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
        
        //if (i==0) {
        //    std::cout << "Position  :  " << pos[0] << " , " << pos[1] << " , " << pos[2] << std::endl;
        //}
    }
   
    n_faces=nf;
    n_tetras=nt;
    
    //number=nv;
    //Export(-2,-2);
    
    return nv;
}
    

// Makes sure everything is in place
void Part_set::GetStarted(){
    //CheckBoxSize();
    //particles.init_neighbour_search(prop->corner_0,prop->corner_1,vbool3(false,false,false),prop->Rsearch);
    //std::cout << "# initiated neighbour serch xith Rsearch" << prop->Rsearch << std::endl;
    
}


// Get the closest neighbours
void Part_set::GetNeighbours() {
   
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



// One simulation step
void Part_set::NextStep(const Simul_props & simul_prop){
    ComputeForces(simul_prop);
    AddConfinementForces(simul_prop);
    IntegrateForces(simul_prop);
    //if (!(prop->elastic)) {
        //particles.update_positions();
    //}
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
void Part_set::Export(int t, double time){
#ifdef HAVE_VTK
    //vtkWriteGrid("particles",t,particles.get_grid(true));
#endif
    ofstream exportfile;
    std::string fname_out="particles_out_"+std::to_string(t)+".txt";
    exportfile.open(fname_out);
    exportfile << "# time : " << time << std::endl;
    exportfile << "# visco : " << prop->visco << std::endl;
    //exportfile << "# pressure : " << simul_prop.pressure << std::endl;
    //myfile << "Writing this to a file.\n";
    exportfile << "# X  Y Z   Nx  Ny  Nz Fx Fy Fz Nneigh ID \n";
    for (int i = 0; i < number; ++i) {
        exportfile <<" "<< get<position>(particles[i])[0] <<" "<< get<position>(particles[i])[1] << " "<< get<position>(particles[i])[2];
        exportfile<< " " << get<orientation>(particles[i])[0] << " " << get<orientation>(particles[i])[1] << " " << get<orientation>(particles[i])[2];
        exportfile<< " " << get<force>(particles[i])[0] << " " << get<force>(particles[i])[1] << " " << get<force>(particles[i])[2];
        exportfile <<" "<< get<nn>(particles[i])<<" " << get<id>(particles[i])<< " \n";
        
    }
    //exampleOutFile.get_comments().push_back("time "+std::to_string(t));
    //exampleOutFile.get_comments().push_back("visco "+std::to_string(prop->visco));
    //exampleOutFile.get_comments().push_back("pressure "+std::to_string(simul_prop.pressure));

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

