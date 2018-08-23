#include "glossary.h"
#include <sstream>
#include "meshless_spin_props.h"


Meshless_props::Meshless_props(const Glossary opt) {
    
    int s = ELASTIC;
    mechanics=s;
    if ( opt.set(s, "mechanics")  ) 
    {
        mechanics=s;
    }
    // Times
    dt=0.01;
    Tend=1.0;
    n_frames=0;
    opt.set(dt,"dt");
    opt.set(Tend,"Tend");
    dt_frames=Tend;
    // Number of times we should save the shape
    opt.set(n_frames,"n_frames");
    nb_times=opt.nb_values("changing_times");
    
    double_list times(nb_times);
    double_list press(nb_times);
    double time=0;
    double pres=0;
    for (int i=0;i<nb_times;++i) {
        opt.set(time ,"changing_times",i);
        opt.set(pres,"changing_pressures",i);
        times[i]=time;
        press[i]=pres;
    }
    
    change_times=times;
    change_press=press;
    
    // Deducing the interval between two saves
    if (n_frames) {
        dt_frames=Tend/(n_frames);
    }
    else {
        dt_frames=Tend;
    }
    n_frames++;
}
