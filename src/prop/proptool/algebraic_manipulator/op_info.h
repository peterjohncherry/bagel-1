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

class Op_Info : public std::enable_shared_from_this<Op_Info>  { 
  friend MultiOp_Info;
 
  private : 
    std::shared_ptr<std::vector<int>> state_ids_;
    char transformation_;

  public :  
    std::string name_;
    std::string op_name_; //TODO temporary hack for CTP names, ctp names should be changed to use state dependent name_

    // should also include symmetry information
     Op_Info( std::string name, std::shared_ptr<std::vector<int>> state_ids, char transformation ) :  
     name_(name), op_name_(name.substr(0,1)),  state_ids_(state_ids), transformation_(tolower(transformation)){}

     Op_Info( std::string name ) : name_(name), op_name_(name.substr(0,1)){}; 
    ~Op_Info(){} 

    virtual std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info() {
      assert( false ); 
      std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> dummy; 
      return dummy;
    }

    virtual std::shared_ptr<Op_Info> op_info( int ii ) { assert( ii == 0 );  return shared_from_this(); }
    virtual int op_order( int ii ) { assert( ii == 0 );  return  0; }
    virtual std::shared_ptr<std::vector<int>> op_order() { std::shared_ptr<std::vector<int>> oo = std::make_shared<std::vector<int>>(1,0); return oo ; }
    virtual char transformation() { return transformation_; } 

};

class MultiOp_Info : public Op_Info { 

  public :  
    int num_ops_;
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> state_ids_;
    std::shared_ptr<std::vector<char>> transformations_;
    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec_;  
    std::shared_ptr<std::vector<int>> op_order_;

    // should also include symmetry information
    MultiOp_Info( std::string name, std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec, const std::vector<int>& op_order ) : 
                   Op_Info( name ), op_info_vec_(op_info_vec), op_order_(std::make_shared<std::vector<int>>(op_order)) {   

      num_ops_ = op_info_vec->size();  
      state_ids_ = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(num_ops_);    
      transformations_ = std::make_shared<std::vector<char>>(num_ops_);    

    } 
    ~MultiOp_Info(){}; 

    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info() { return op_info_vec_; }
    std::shared_ptr<Op_Info> op_info( int ii ) { return (*op_info_vec_)[ii]; }

    std::shared_ptr<std::vector<int>> op_order() { return op_order_; }
    int op_order( int ii ) { return (*op_order_)[ii]; }

    char transformation() { throw std::logic_error("Should not try to transform MultiOp all as one! Aborting!! "); return 'X'; } 
};

#endif


//     std::vector<std::vector<int>>::iterator si_it = state_ids_->begin();
//     std::vector<char>::iterator t_it =  transformations_->begin();
//     for ( std::vector<std::shared_ptr<Op_Info>>::iterator oi_it = op_info_vec->begin(); oi_it != op_info_vec->end(); oi_it++, si_it++, t_it++ ) { 
//        *t_it = (*oi_it)->transformation_;
//        *si_it =  (*oi_it)->state_ids_;
//      }

//    std::shared_ptr<std::vector<int>> state_ids(int ii){ assert( ii == 0 ); return (*state_ids_)[ii]; };
//    char transformations( int ii ) { assert( ii == 0 ); return transformation;  };
