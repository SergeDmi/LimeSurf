/*
  This is part of Limesurf
    A program meant to integrated the shape of triangulated surface under constraints.
    
    Copyright 2017-2020 Serge Dmitrieff, Institut Jacques Monod, CNRS
    www.biophysics.fr
*/
#include "mesh.h"
#include "yaml-cpp/yaml.h" 
using namespace std;

/*
 Definitions 
*/

// A vector of simulation properties
typedef std::vector<Simul_props> properties;
// A vector of meshes
typedef std::vector<Mesh>        meshes    ;

/*
 Help functions
*/

// Builds the vector of properties from the config file
void const make_props(properties & props, const YAML::Node & yaconf) {
    YAML::Node all_runs = yaconf["runs"];
    Simul_props prop;
    for(YAML::const_iterator it=all_runs.begin();it!=all_runs.end();++it) {
        // We keep reading in the same property
        // This means that we keep by default the previous simul properties, unless overwritten by Read_config
        prop.Read_config(it);
        props.push_back(prop);
    }
}

// Builds the vector of meshes from the config file
void const make_meshes(meshes & meshuggas, const YAML::Node & yaconf) {
    YAML::Node all_meshes = yaconf["meshes"];
    Part_set_props prop;
    for(YAML::const_iterator it=all_meshes.begin();it!=all_meshes.end();++it) {
        // We keep reading in the same property
        // This means that we keep by default the previous wall properties, unless overwritten by Read_config
        // Warning : only generic properties ( Part_set_props) are kept
        // Specialized properties ( e.g. Elastic_set_props ) are not kept
        prop.Read_config(it);
        Mesh mesh(it,prop);
        // Mesh initiation 
        mesh.Initiate();
        // Push back into meshuggas. Djent. 
        meshuggas.push_back(mesh);
    }
}

/*
 Part where we actually run the simulations
*/

// Make_run runs a simulation of properties prop for all meshes in meshuggas. 
int make_run(const Simul_props & prop, meshes & meshuggas, double t0) {
    // Program state
    int dvg=0;
    //std::cout << "Current pressure " << prop.pressure << std::endl;
    // Running all the meshes	
    for (auto &mesh: meshuggas) {
        // Preparing times
        double t_save=0;
        int n_save=0;
        double t=0;
        bool is_diverging=0;
        mesh.mesh_parts->Summary();
        // Runnning until Tend, or divergence
        while (t<prop.Tend && !is_diverging) {
            t+=prop.dt;
            // the mesh performs its next stp
            mesh.mesh_parts->NextStep(prop);
            is_diverging=mesh.mesh_parts->is_diverging();
            // Saving if needed
            t_save+=prop.dt;
            if (t_save>prop.dt_frames ) {
                t_save=0;
                mesh.mesh_parts->Export_bly(n_save,prop,t+t0);
                n_save++;
                mesh.mesh_parts->Summary();
            }
        }
        // implicit bool to int conversion for the win (and less branching)
        mesh.mesh_parts->Summary();
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
   
    // No export here but convenient to keep the code ;)
    //all_meshes[0].cell_wall->Export(-1,-1);
   
    // For all simulation runs (i.e. all sets of properties) we simulate
	double time=0.0;
    for (const auto &run: runs) {
        dvg=make_run(run,all_meshes,time);
		time+=run.Tend;
     }
	 
    // Returning program state
    return dvg;
}


