#include "spin_parts.h"
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
    
    //glos.read_strings(argc-1, argv+1);
    std::cout << "# Created set1" << std::endl;
    if (argc>2) {
        std::cout << "# +++Created set1" << std::endl;
        std::string fname(argv[2]);
        set1.create(fname);
        std::cout << "# Created set1--++++++" << std::endl;
    }
    else
    {
        std::cout << "# -- Created set1" << std::endl;
        set1.create();
        std::cout << "# Created set1-----" << std::endl;
    }


    set1.GetStarted();
    double t=0;
    while (t<simul_prop.Tend) {
        std::cout << "#"  << std::flush;
        t+=simul_prop.dt;
        set1.NextStep(&simul_prop);
        
    }
    std::cout << std::endl;
    set1.ComputeForces();
    //set1.GetNeighbours();
    set1.Export(0);
    return N;
}
