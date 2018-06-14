//#include "spin_parts.h"
#include "elastic_spin_parts.h"
#include <sstream>
#include "glossary.h"
#include <fstream>


int main(int argc, char* argv[])
{
    //std::istringstream fname(argv[1]);
    ifstream config(argv[1]);
    Glossary glos(config);
    Meshless_props simul_prop(glos);
    Part_set_props set_prop(glos);
    Part_set set1(&set_prop);
    Elastic_part_set set2(&set_prop);
    std::string fname_in="cell.ply";
    std::string fname_out="simulated_cell.ply";
    //glos.read_strings(argc-1, argv+1);
    std::cout << "# Created set1" << std::endl;
    if (argc>2) {
        //std::cout << "# +++Created set1" << std::endl;
        fname_in=argv[2];
        set2.create(fname_in);
        //set1.GetStarted();
        //set1.Export_bly(fname);
        //std::cout << "# Created set1--++++++" << std::endl;
    }
    else
    {
        //std::cout << "# -- Created set1" << std::endl;
        set2.create();
        //std::cout << "# Created set1-----" << std::endl;
    }


    set2.GetStarted();
    double t=0;
    while (t<simul_prop.Tend && !set2.diverging) {
        //std::cout << "#"  << std::flush;
        t+=simul_prop.dt;
        set2.NextStep(&simul_prop);
        
    }
    std::cout << std::endl;
    set2.ComputeForces();
    set2.Export(0);
    set2.Export_bly(fname_out);
    return (int)set2.diverging;
}
