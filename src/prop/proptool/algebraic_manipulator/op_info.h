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

class MultiOp_Info;

class Op_Info { 
  friend MultiOp_Info;
 
  private : 
    std::shared_ptr<std::vector<int>> state_ids_;
    std::shared_ptr<std::vector<int>> state_ids(){ return state_ids_; }
    char transformation_;

  public :  
    std::string name_;

    // should also include symmetry information
     Op_Info( std::string name, std::shared_ptr<std::vector<int>> state_ids, char transformation ) :  
     name_(name), state_ids_(state_ids), transformation_(transformation) {}
     Op_Info( std::string name ) : name_(name) {}; 
    ~Op_Info(){} 

};

class MultiOp_Info : public Op_Info { 

  public :  
    int num_ops_;
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> state_ids_;
    std::shared_ptr<std::vector<char>> transformations_;
    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec_;  

    // should also include symmetry information
    MultiOp_Info( std::string name, std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec ) : 
                   Op_Info( name ) {   

      op_info_vec_ =  op_info_vec;
      num_ops_ = op_info_vec->size();  
      state_ids_ = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(num_ops_);    
      transformations_ = std::make_shared<std::vector<char>>(num_ops_);    

//      std::vector<std::vector<int>>::iterator si_it = state_ids_->begin();
//      std::vector<char>::iterator t_it =  transformations_->begin();
//     for ( std::vector<std::shared_ptr<Op_Info>>::iterator oi_it = op_info_vec->begin(); oi_it != op_info_vec->end(); oi_it++, si_it++, t_it++ ) { 
//        *t_it = (*oi_it)->transformation_;
//        *si_it =  (*oi_it)->state_ids_;
//      }

    } 

   ~MultiOp_Info(){}; 

};


#endif
