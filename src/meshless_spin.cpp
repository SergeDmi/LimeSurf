//#include "spin_parts.h"
#include "simple_viscoel_spin_parts.h"
#include "simple_viscoel_parts_props.h"
#include "meshless_spin_props.h"
#include <sstream>
#include "glossary.h"
#include <fstream>


int main(int argc, char* argv[])
{
    // Reading input
    ifstream config(argv[1]);
    
    // Making a glossary from config file
    Glossary glos(config);
    
    // simulation & properties 
    Meshless_props simul_prop(glos);
    
    // Empty creation of particle set and properties
    Part_set * cell_wall;
    
    // Now we properly define our particle set
    switch (simul_prop.mechanics) {
        case 0 : {
            Part_set_props wall_props(glos);
            Part_set * wall=new Part_set(&wall_props);
            cell_wall=wall;
            std::cout << "# Starting viscous simulation" << std::endl;
            break; }
        case 1 : {  
            Elastic_set_props wall_props(glos);
            Elastic_part_set *wall = new Elastic_part_set(&wall_props);
            cell_wall=wall;
            std::cout << "# Starting elastic simulation" << std::endl;
            break; }
        case 2 : {
            Simple_viscoel_set_props wall_props(glos);
            Simple_viscoel_part_set * wall = new Simple_viscoel_part_set(&wall_props);
            cell_wall=wall;
            std::cout << "# Starting viscoelastic simulation" << std::endl;
            break; }
    }
    // Properties of the particle set
    //Elastic_set_props set_prop(glos);
    
    //Part_set set1(&set_prop);
    // Now making an elastic simulation
    // @TODO : this should be done much more elegantly
    //Elastic_part_set set2(&set_prop);
    //Part_set test_set(&set_prop);
    //test_set=set2; 
   
    // Export name
    
  
    cell_wall->create();   
    

    std::string fname_out="simulated_cell";
    // Getting ready to simulate
    cell_wall->GetStarted();
    
    // Preparing times
    double t_save=0;
    int n_save=0;
    double t=0;
    int n_change=0;
    int n_export=0;
    
    // We export simulation state in a text file
    cell_wall->Export(n_export);
    n_export++;
    cell_wall->Export_bly(fname_out,n_save,&simul_prop);
    n_save++;
    
    // Simulation steps
    while (t<simul_prop.Tend && !cell_wall->is_diverging()) {
        t+=simul_prop.dt;
        cell_wall->NextStep(&simul_prop);
        
        // Saving if needed
        t_save+=simul_prop.dt;
        if (t_save>simul_prop.dt_frames ) {
            t_save=0;
            cell_wall->Export_bly(fname_out,n_save,&simul_prop);
            n_save++;
            
        }
        
       // Changing pressure if specified in the config file
       // we export sim state at each pressure change
       if (n_change<simul_prop.nb_times) {
           if (t>simul_prop.change_times[n_change]) {
               cell_wall->prop->pressure=simul_prop.change_press[n_change];
               std::cout << "# Changing pressure to " << cell_wall->prop->pressure << std::endl;
               cell_wall->Export(n_export);
               n_export++;
               n_change++;
           }
       }
    }
    
    // Wrapping up
    std::cout << std::endl;
    cell_wall->ComputeForces();
    cell_wall->Export(n_export);
    cell_wall->Export_bly(fname_out,n_save,&simul_prop);
    return (int)cell_wall->is_diverging();
}
