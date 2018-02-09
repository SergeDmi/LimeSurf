#include <random>
#include "glossary.h"
#include "Aboria.h"
using namespace Aboria;
#include <sstream>
#include "spin_parts_props.h"


Part_set_props::Part_set_props(const Glossary opt) {
    init_shape=0;
    init_number=100;
    init_radius=50;
    Rmax=3;
    minR=0.5;
    L=200.0;
    k_align=1.0;
    k_att=2.0;
    k_rep=1.0;
    k_bend=1.0;
    p_rep=12.0;
    p_att=6.0;
    visco=1.0;
    Rvisc=1.0;
    opt.set(init_shape, "shape", KeyList<int>("sphere", 0, "pombe", 1, "sheet", 2));
    opt.set(init_number, "number");
    opt.set(init_radius, "radius");
    opt.set(visco, "visco");
        opt.set(Rvisc, "Rvisco");
    opt.set(Rmax, "Rmax");
    opt.set(minR, "Rmin");
    opt.set(L, "box");
    R0=0.1;
}
