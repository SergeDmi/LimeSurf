#include <random>
#include "Aboria.h"
using namespace Aboria;
#include "spin_parts_props.h"
#include "meshless_spin_props.h"
#include <boost/math/constants/constants.hpp>
#include <math.h>
#define MAXNEIGHBOURS 30


class Part_set {
    
    std::vector<int> NEIGHBOURS;
    ABORIA_VARIABLE(orientation,vdouble3,"orientation")
    ABORIA_VARIABLE(torque,vdouble3,"orientation")
    ABORIA_VARIABLE(force ,vdouble3,"orientation")
    ABORIA_VARIABLE(neighbours,std::vector<int>,"neighbours")
    ABORIA_VARIABLE(nn,double,"neighbour number")
    ABORIA_VARIABLE(state,double,"state")
    
    typedef Particles<std::tuple<orientation,neighbours,force,torque,nn,state>,3> particle_type;
    //typedef Particles<std::tuple<orientation>,<test>,2> container_type;
    
    typedef typename particle_type::position position;
public:
    Part_set(Part_set_props*);
    void create();
    void create(std::string);
    Part_set_props* prop;
    int num();
    //int PutOnSphere(int,double);
    void NextStep(const Meshless_props*);
    void IntegrateForces(const Meshless_props*);
    void ComputeForcesViscous();
    void ComputeForcesElastic();
    void ComputeForces();
    void GetNeighbours();
    void Export(int);
    void GetStarted();
    void ClearForces();
    int PutOnSphere();
    int PutOnSheet();
    int load_from_text(std::string);
    int max_neighbours;
private:
    particle_type particles;
    int number;

};
