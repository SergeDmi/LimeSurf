#include "glossary.h"
#include <sstream>
#include "meshless_spin_props.h"


Meshless_props::Meshless_props(const Glossary opt) {
    // Times
    dt=0.01;
    Tend=1.0;
    n_frames=0;
    opt.set(dt,"dt");
    opt.set(Tend,"Tend");
    dt_frames=Tend;
    // Number of times we should save the shape
    opt.set(n_frames,"n_frames");
    
    // Deducing the interval between two saves
    if (n_frames) {
        dt_frames=Tend/n_frames;
    }
    else {
        dt_frames=Tend;
    }
    n_frames++;
}
