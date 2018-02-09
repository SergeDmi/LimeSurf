#include <random>
#include "glossary.h"
#include "Aboria.h"
using namespace Aboria;
#include <sstream>
#include "spin_parts_props.h"


Part_set_props::Part_set_props(const Glossary opt) {
    init_shape=1;
    init_number=1;

    Rmax=0.3;
    minR=0.05;
    L=4.0;
     k_align=1.0;
    k_att=1.0;
        k_rep=1.0;
    k_bend=1.0;
    opt.set(init_shape, "shape", KeyList<int>("sphere", 0, "pombe", 1, "sheet", 2));
    opt.set(init_number, "number");
    opt.set(R0, "R0");
    opt.set(k_att, "attraction");
    R0=0.1;
}
