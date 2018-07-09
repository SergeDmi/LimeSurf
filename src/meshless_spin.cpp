//#include "spin_parts.h"
#include "elastic_spin_parts.h"
#include "elastic_parts_props.h"
#include <sstream>
#include "glossary.h"
#include <fstream>


int main(int argc, char* argv[])
{
    // Reading input
    ifstream config(argv[1]);
    
    // Making a glossary from config file
    Glossary glos(config);
    
    // simulation properties
    Meshless_props simul_prop(glos);
    
    // Properties of the particle set
    Elastic_set_props set_prop(glos);
    
    //Part_set set1(&set_prop);
    // Now making an elastic simulation
    // @TODO : this should be done much more elegantly
    Elastic_part_set set2(&set_prop);
    
    // Dummy inpyt name
    std::string fname_in="cell.ply";
    
    // Export name
    std::string fname_out="simulated_cell";
    
    if (argc>2) {
        
        fname_in=argv[2];
        set2.create(fname_in);
        
    }
    else
    {
        
        set2.create();
        
    }

    
    // Getting ready to simulate
    set2.GetStarted();
    
    // Preparing times
    double t_save=0;
    int n_save=0;
    double t=0;
    
    // Simulation steps
    while (t<simul_prop.Tend && !set2.is_diverging()) {
        t+=simul_prop.dt;
        set2.NextStep(&simul_prop);
        
        // Saving if needed
        t_save+=simul_prop.dt;
        if (t_save>simul_prop.dt_frames) {
            t_save=0;
            set2.Export_bly(fname_out,n_save,&simul_prop);
            n_save++;
            
        }
        
    }
    std::cout << std::endl;
    set2.ComputeForces();
    set2.Export(0);
    set2.Export_bly(fname_out,n_save,&simul_prop);
    return (int)set2.is_diverging();
}
