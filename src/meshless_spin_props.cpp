#include "glossary.h"
#include <sstream>
#include "meshless_spin_props.h"


Meshless_props::Meshless_props(const Glossary opt) {
    dt=0.1;
    opt.set(dt,"dt");

}
