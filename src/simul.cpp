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
#include "glossary.h"
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
    Mesh mesh;
    for(YAML::const_iterator it=all_meshes.begin();it!=all_meshes.end();++it) {
        mesh.Read_config(it);
        meshugas.push_back(mesh);
    }
}

void make_run(const Simul_props & prop, meshes meshugas) {
    // Preparing times
    double t_save=0;
    int n_save=0;
    double t=0;
    int n_change=0;
    int n_export=0;
    bool is_diverging=0;
    for (auto mesh: meshugas) {
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
    }
    
    
}

int main(int argc, char* argv[])
{
    // Reading input
    ifstream config(argv[1]);
    // Making a glossary from config file
    Glossary glos(config);
    //YAML::Node yaconf;
    //YAML::Node yaconf2 = YAML::LoadFile("conf.yaml");
    YAML::Node yaconf2 = YAML::LoadFile(argv[2]);
    Mesh mesh(glos);
    meshes all_meshes;
    mesh.Initiate();
    //Simul_props simul_prop(glos);
    properties runs;
    make_props(runs,yaconf2);
    //make_meshes(all_meshes,yaconf2);
    //mesh
    all_meshes.push_back(mesh);
    std::cout << "fname in should in meshless : " << mesh.cell_wall->prop->fname_in << std::endl;
    mesh.cell_wall->create();   
    std::cout << "Created cell wall " << std::endl;
    
    std::string fname_out="simulated_cell";
    // Getting ready to simulate
    mesh.cell_wall->GetStarted();
    
    
    
    //for(auto run=runs.begin();run!=runs.end();++run) {
    for (auto run: runs) {
        make_run(run,all_meshes);
     }
     
    /*
    Simul_props simul_prop=runs[0];
    
    YAML::Node buns = yaconf2["runs"];
    int a=0;
     */
     /*
    for(auto it=buns.begin();it!=buns.end();++it) {
    
        std::cout << "Run # " <<  a << "as" << it->first.as<std::int8_t>() <<  std::endl;
        YAML::Node conf=it->second;
        if (conf["pressure"])
            std::cout << "         pressure : " << conf["pressure"].as<std::string>() <<  std::endl;
        if (conf["test"])
            std::cout << "         test : " << conf["test"].as<std::string>() <<  std::endl;
        //std::cout << "         pressure : " << conf["pressure"].as<std::string>() <<  std::endl;
        //for(auto jt=it->second.begin();jt!=it->second.end();++jt) {
        //    std::cout << "??? " << jt->first.as<std::string>() <<  std::endl;
        //}
        a++;
        //mesh.Update(*it);
        
    }
     */
    
    
    // simulation & properties 
    
    /*
    // Empty creation of particle set and properties
    Part_set * cell_wall;
    Part_set_props * wall_props;
    
    // Now we properly define our particle set
    switch (simul_prop.mechanics) {
        case 0 : {
            Part_set_props  * props = new Part_set_props(glos);
            wall_props=props;
            Part_set * wall=new Part_set(props);
            cell_wall=wall;
            std::cout << "# Starting viscous simulation" << std::endl;
            break; }
        case 1 : {  
            Elastic_set_props  * props = new Elastic_set_props(glos);
            wall_props=props;
            Elastic_part_set * wall=new Elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting elastic simulation" << std::endl;
            //std::cout << "fname in should in meshless : " << wall_props->fname_in << std::endl;
            break; }
        case 2 : {
            Simple_viscoel_set_props *props = new Simple_viscoel_set_props(glos);
            wall_props=props;
            Simple_viscoel_part_set * wall = new Simple_viscoel_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting viscoelastic simulation" << std::endl;
            break; }
        case 3 : {
            Triangle_set_props  * props = new Triangle_set_props(glos);
            wall_props=props;
            Triangle_part_set * wall=new Triangle_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting triangular elastic simulation" << std::endl;
            break; }                           
        case 4 : {
            Elastic_set_props  * props = new Elastic_set_props(glos);
            wall_props=props;
            Tetr_elastic_part_set *wall = new Tetr_elastic_part_set(props);
            cell_wall=wall;
            std::cout << "# Starting 3D elastic simulation" << std::endl;
            break; }
    }
    
    
    */
    // Properties of the particle set
    //Elastic_set_props set_prop(glos);
    
    //Part_set set1(&set_prop);
    // Now making an elastic simulation
    // @TODO : this should be done much more elegantly
    //Elastic_part_set set2(&set_prop);
    //Part_set test_set(&set_prop);
    //test_set=set2; 
   
    // Export name
    
    
    
    
    // We export simulation state in a text file
    //mesh.cell_wall->Export(n_export);
    //n_export++;
    //mesh.cell_wall->Export_bly(fname_out,n_save,&simul_prop);
    //n_save++;
    
    //std::cout << "Starting sim " << std::endl;
    // Simulation steps
    /*
    while (t<simul_prop.Tend && !mesh.cell_wall->is_diverging()) {
        t+=simul_prop.dt;
        //std::cout << t << std::endl;
        mesh.cell_wall->NextStep(&simul_prop);
        
        // Saving if needed
        t_save+=simul_prop.dt;
        if (t_save>simul_prop.dt_frames ) {
            t_save=0;
            mesh.cell_wall->Export_bly(fname_out,n_save,&simul_prop);
            n_save++;
            
        }
        */
       // Changing pressure if specified in the config file
       // we export sim state at each pressure change
       /*
       if (n_change<simul_prop.nb_times) {
           if (t>simul_prop.change_times[n_change]) {
               mesh.cell_wall->prop->pressure=simul_prop.change_press[n_change];
               std::cout << "# Changing pressure to " << mesh.cell_wall->prop->pressure << std::endl;
               mesh.cell_wall->Export(n_export);
               n_export++;
               n_change++;
           }
       }
        */
    
    
    // Wrapping up
    //std::cout << std::endl;
    /*
    mesh.cell_wall->ComputeForces(simul_prop);
    mesh.cell_wall->Export(n_export);
    mesh.cell_wall->Export_bly(fname_out,n_save,&simul_prop);
     */
    mesh.cell_wall->Export(0);
    return (int)mesh.cell_wall->is_diverging();
}
