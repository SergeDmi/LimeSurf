#include <random>
#include "Aboria.h"
using namespace Aboria;
#include "spin_parts_props.h"
#include "meshless_spin_props.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>
#define MAXNEIGHBOURS 100


class Part_set {
    
    typedef Vector<int,100> NEIGHBOURS;
    ABORIA_VARIABLE(orientation,vdouble3,"orientation")
    ABORIA_VARIABLE(neighbours,NEIGHBOURS,"neighbours")
    ABORIA_VARIABLE(nn,double,"neighbour number")
    ABORIA_VARIABLE(state,double,"state")
    
    typedef Particles<std::tuple<orientation,neighbours,nn,state>,3> particle_type;
    //typedef Particles<std::tuple<orientation>,<test>,2> container_type;
    
    typedef typename particle_type::position position;
public:

    Part_set(Part_set_props*);
    Part_set_props* prop;
    int num();
    int PutOnSphere(int,double);
    void NextStep(Meshless_props*);
    void Export(int);
    void GetStarted();
private:
    particle_type particles;
    int number;
};