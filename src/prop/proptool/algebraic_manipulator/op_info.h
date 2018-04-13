#ifndef __SRC_PROP_PROPTOOL_Op_Info_H
#define __SRC_PROP_PROPTOOL_Op_Info_H

#include <vector>
#include <memory>
#include <string>

// Small class to label state specfic operator and connected symmetries
// Necessary for generation of appropriate range blocks.
// Information for specification is output during initialization of equation.
// This, or something very similar should be the key to the MT_map.
// Likewise for CTP; we need state information, which is currently missing.

class Op_Info { 
  public :  

    std::vector<int> state_ids_;
    std::string name_;
    // should also include symmetry information
     Op_Info( std::string name, const std::vector<int>& state_ids ) :  name_(name), state_ids_(state_ids) {}
    ~Op_Info(){} 

};

#endif
