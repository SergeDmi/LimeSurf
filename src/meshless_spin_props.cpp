#include "glossary.h"
#include <sstream>
#include "meshless_spin_props.h"


Meshless_props::Meshless_props(const Glossary opt) {
    dt=0.01;
    Tend=1.0;
    opt.set(dt,"dt");
    opt.set(Tend,"Tend");
    opt.set(n_frames,"n_frames");
    if (n_frames) {
        dt_frames=Tend/n_frames;
    }
    else {
        dt_frames=Tend;
    }
    n_frames++;
}
