/*
 #include "spin_parts.h"
#include "tetr_elastic_spin_parts.h"
#include "simple_viscoel_spin_parts.h"
#include "simple_viscoel_parts_props.h"
#include "meshless_spin_props.h"
#include "elastic_triangle_parts.h"
#include "triangle_parts_props.h"
 */
#include "mesh.h"
#include <sstream>
//#include "glossary.h"
#include <fstream>
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace std;

typedef std::vector<Simul_props> properties;
typedef std::vector<Mesh>        meshes    ;


void const make_props(properties & props, const YAML::Node & yaconf) {
    YAML::Node all_runs = yaconf["runs"];
    Simul_props prop;
    for(YAML::const_iterator it=all_runs.begin();it!=all_runs.end();++it) {
        prop.Read_config(it);
        props.push_back(prop);
    }
}

void const make_meshes(meshes & meshugas, const YAML::Node & yaconf) {
    YAML::Node all_meshes = yaconf["meshes"];
   
    int a=0;
    for(YAML::const_iterator it=all_meshes.begin();it!=all_meshes.end();++it) {
        //std::cout << "Making mesh # " << a << std::endl;
        Mesh mesh(it);
        mesh.Initiate();
        //mesh.cell_wall->create();   
        //mesh.Read_config(it);
        meshugas.push_back(mesh);
        //std::cout << "Made mesh ! " << std::endl;
        a++;
        //mesh.
        mesh.cell_wall->Export(-2,-2);
    }
}

int make_run(const Simul_props & prop, meshes & meshugas) {
    // Preparing times
    int dvg=0;
    for (auto mesh: meshugas) {
        double t_save=0;
        int n_save=0;
        double t=0;
        int n_change=0;
        int n_export=0;
        bool is_diverging=0;
        //mesh.cell_wall->GetStarted();
        mesh.cell_wall->Export(0,0);
        while (t<prop.Tend && !is_diverging) {
            t+=prop.dt;
            //std::cout << t << std::endl;
            mesh.cell_wall->NextStep(prop);
            is_diverging=mesh.cell_wall->is_diverging();
            // Saving if needed
            t_save+=prop.dt;
            if (t_save>prop.dt_frames ) {
                t_save=0;
                mesh.cell_wall->Export_bly(n_save,prop,t);
                n_save++;
                
            }
        }
        if (is_diverging) {
            dvg++;
        }
        mesh.cell_wall->Export(1,t);
    }
    
    return dvg;
}

int main(int argc, char* argv[])
{
    // Reading input
    int dvg=0;
    //ifstream config(argv[1]);
    // Making a glossary from config file
    //Glossary glos(config);
    YAML::Node yaconf2 = YAML::LoadFile(argv[1]);
    
    
    
    
    
    // Making all properties
    properties runs;
    make_props(runs,yaconf2);
    
    // Making all meshes
    meshes all_meshes;
    make_meshes(all_meshes,yaconf2);
   
    all_meshes[0].cell_wall->Export(-1,-1);
   
    // Running and gunning
    for (auto run: runs) {
        dvg+=make_run(run,all_meshes);
        all_meshes[0].cell_wall->Export(-3,-3);
     }
     
    return dvg;
}
