#include <random>
#include "glossary.h"
//#include "Aboria.h"
//using namespace Aboria;
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
    R0=1.0;
    e=1.0;
    Rvisc=1.0;
    elastic=0;
    Fmax=DBL_MAX;
    pressure=0;
    max_neighbours=9;
    double box[6];
    
    opt.set(init_shape, "shape", KeyList<int>("sphere", 0, "sheet", 1, "pombe", 2));
    opt.set(init_number, "number");
    opt.set(init_radius, "radius");
    
    opt.set(visco, "visco",0);
    opt.set(Rvisc, "visco",1);
    opt.set(Rvisc, "Rvisco");
    opt.set(max_neighbours,"max_neighbours");
    
   
    opt.set(k_att, "attraction",0);
    opt.set(k_rep, "repulsion",0);
    opt.set(p_att, "attraction",1);
    opt.set(p_rep, "repulsion",1);
    
    p_bend=p_att;
    p_align=p_att;
    
    opt.set(k_align, "align",0);
    opt.set(k_bend, "bending",0);
    opt.set(p_align, "align",1);
    opt.set(p_bend, "bending",1);
    opt.set(e, "prestrain");
    opt.set(elastic, "elastic");
    
    for (int i=0; i<6; ++i) {
        opt.set(box[i], "box",i);
    }
    

    corner_0=vdouble3(box[0],box[1],box[2]);
    corner_1=vdouble3(box[3],box[4],box[5]);
    std::cout << "# bounding box from " << corner_0 << " to " << corner_1 << std::endl;
    std::cout << "# prestrain " << e << std::endl;
    // Computing R0
    if (k_att>0) {
            R0=pow(p_rep*k_rep/(p_att*k_att),1.0/(p_rep-p_att));
    }
    opt.set(R0, "R0");
    minR=0.75*R0;
    Rmax=2.0*R0;

    k_elast=1.0;
    opt.set(minR, "Rmin");
    minR2=pow(minR,2.0);
    if (k_rep>0) {
        Fmax=k_rep/pow(minR,p_rep);
    }
    opt.set(Fmax, "Fmax");
    Fmax2=Fmax*Fmax;
    //opt.set(L, "box");
    opt.set(Rmax, "Rmax");
    Rsearch=Rmax;
    opt.set(Rsearch, "Rsearch");
    opt.set(Rmax, "Rmax");
    
    opt.set(pressure, "pressure");
    // Adimentionalized pressure
    pressure*=sqrt(3.0)/(2.0);
    
    std::cout << "# Fmax2=" << Fmax2 << std::endl;
}
