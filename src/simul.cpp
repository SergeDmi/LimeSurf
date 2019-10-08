#include "mesh.h"
#include "yaml-cpp/yaml.h"  // IWYU pragma: keep

using namespace std;
// A vector of simulation properties
typedef std::vector<Simul_props> properties;
// A vector of meshes
typedef std::vector<Mesh>        meshes    ;


// Builds the vector of properties from the config file
void const make_props(properties & props, const YAML::Node & yaconf) {
    YAML::Node all_runs = yaconf["runs"];
    Simul_props prop;
    for(YAML::const_iterator it=all_runs.begin();it!=all_runs.end();++it) {
        prop.Read_config(it);
        props.push_back(prop);
    }
}

// Builds the vector of meshes from the config file
void const make_meshes(meshes & meshugas, const YAML::Node & yaconf) {
    YAML::Node all_meshes = yaconf["meshes"];
    for(YAML::const_iterator it=all_meshes.begin();it!=all_meshes.end();++it) {
        Mesh mesh(it);
        mesh.Initiate();
        meshugas.push_back(mesh);
    }
}

int make_run(const Simul_props & prop, meshes & meshugas) {
    // Program state
    int dvg=0;
    // Running all the meshes
    for (auto mesh: meshugas) {
        // Preparing times
        double t_save=0;
        int n_save=0;
        double t=0;
        bool is_diverging=0;
        // Runnning for all times 
        while (t<prop.Tend && !is_diverging) {
            t+=prop.dt;
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
        // implicit bool to int conversion for the win (and less branching)
        dvg+=is_diverging;
    }
    
    return dvg;
}



int main(int argc, char* argv[])
{
    // Program state
    int dvg=0;
    
    // Reading input
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
     }
     
    // Returning program stage
    return dvg;
}


