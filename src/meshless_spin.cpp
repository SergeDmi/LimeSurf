#include "spin_parts.h"
#include <sstream>
#include "glossary.h"
#include <fstream>


int main(int argc, char* argv[])
{
    //std::istringstream fname(argv[1]);
    ifstream config(argv[1]);
    Glossary glos(config);
    //glos.read_strings(argc-1, argv+1);
    Part_set_props set_prop(glos);
    Meshless_props simul_prop(glos);

    Part_set set1(&set_prop);
    
    set1.GetStarted();
    set1.NextStep(&simul_prop);
    set1.Export(0);
    return N;
}