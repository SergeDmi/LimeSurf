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
    opt.set(init_shape, "shape", KeyList<int>("sphere", 0, "sheet", 1, "pombe", 2));
    opt.set(init_number, "number");
    opt.set(init_radius, "radius");
    opt.set(visco, "visco");
    opt.set(Rvisc, "Rvisco");
    opt.set(Rmax, "Rmax");
    opt.set(k_align, "align");
    opt.set(k_bend, "rigidity");
    
    R0=pow(p_rep*k_rep/(p_att*k_att),1.0/(p_rep-p_att));
    minR=0.75*R0;
    minR2=pow(minR,2.0);
    Fmax=k_rep*pow(minR,p_rep);
    opt.set(minR, "Rmin");
    opt.set(L, "box");
    
}
