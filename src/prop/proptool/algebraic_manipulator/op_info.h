#ifndef __SRC_PROP_PROPTOOL_Op_Info_H
#define __SRC_PROP_PROPTOOL_Op_Info_H

#include <vector>
#include <memory>
#include <string>

#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>

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
    std::string op_name_;
    std::string op_state_name_;
    std::string op_full_name_;
    std::string name_;
 
    std::shared_ptr< std::vector<std::pair<std::vector<int>,std::shared_ptr<Transformation>>> > equivalent_operators_list_; 

    bool spin_symm_;
    std::pair<double,double> spin_flip_factor_ ;

    // TODO Find way round putting these in base, preferably no involving using functions
    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec_;  
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> state_ids_list_;
    std::shared_ptr<std::vector<char>> transformations_;
    std::shared_ptr<std::vector<int>> op_order_;


    // should also include symmetry information
     Op_Info( std::string op_name,  std::string op_state_name, std::string op_full_name, std::shared_ptr<std::vector<int>> state_ids, char transformation ) :  
     op_name_(op_name), op_state_name_(op_state_name), op_full_name_(op_full_name), name_(op_full_name), state_ids_(state_ids), transformation_(tolower(transformation)){

       std::cout << "op_state_name = " << op_state_name_ << std::endl;
     }

     Op_Info( std::string op_name, std::string op_state_name, std::string op_full_name ,
              std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec, const std::vector<int>& op_order ) :
              op_name_(op_name), op_state_name_( op_state_name ), op_full_name_( op_full_name ), name_( op_full_name ),
              op_info_vec_(op_info_vec), op_order_(std::make_shared<std::vector<int>> (op_order) ){
              std::cout << "op_state_name = " << op_state_name_ << std::endl;
    }; 
  
    ~Op_Info(){} 

    //TODO This is awkward, should probably always have these, but how to get round inability to set shared from this in constructor?
    //     Don't want to use inline function every time I access op_info_vec/op_order/transformations etc.
    void initialize_as_multiop_info(){
        op_info_vec_ = std::make_shared<std::vector<std::shared_ptr<Op_Info>>>( 1, shared_from_this()); 
        transformations_ = std::make_shared<std::vector<char>>( 1, transformation_ ); 
        op_order_ = std::make_shared<std::vector<int>>( 1, 0 ); 
        state_ids_list_ = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>> (1, state_ids_);
        return;
    } ; 

    virtual std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info() {
      assert( false ); 
      std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> dummy; 
      return dummy;
    }

    std::shared_ptr<std::vector<char>> transformations() { return transformations_; } 
    std::shared_ptr<std::vector<int>> op_order() { op_order_; }

    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec() { return op_info_vec_; }
    std::shared_ptr<Op_Info> op_info( int ii ) { return (*op_info_vec_)[ii]; }
    std::shared_ptr<Op_Info> op_info_canonical( int ii ) { return (*op_info_vec_)[(*op_order_)[ii]]; } 

    virtual int op_order( int ii ) { assert( ii == 0 );  return  0; }
    virtual char transformation() { return transformation_; } 
    virtual std::string op_state_name_canonical() { return op_state_name_; } 

};

//Note : The order of members in op_info_vec is determined by the corresponding braket.
//       Functions with "canonical" will return names/op_info as though they were in the "canonical" order; that is the
//       same order which all MultiTensOps, and all data structures used in tensor arithmetic, are built.
class MultiOp_Info : public Op_Info { 

  public :  
    int num_ops_;
    std::string op_state_name_canonical_;

    // TODO currently fudged; would like to only have op_info_vec defined in MultiOp_Info class
    MultiOp_Info( std::string op_name, std::string op_state_name, std::string op_full_name, 
                  std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec, const std::vector<int>& op_order ) : 
                  Op_Info( op_name, op_state_name, op_full_name, op_info_vec, op_order ) {

      num_ops_ = op_info_vec->size();  
      state_ids_list_ = std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>(num_ops_);    

      transformations_ = std::make_shared<std::vector<char>>(num_ops_);    
      std::vector<char>::iterator t_it = transformations_->begin(); 
      std::vector<std::shared_ptr<std::vector<int>>>::iterator sil_it = state_ids_list_->begin(); 
      for ( std::vector<std::shared_ptr<Op_Info>>::iterator oi_it = op_info_vec_->begin(); oi_it != op_info_vec_->end(); oi_it++ , t_it++, sil_it++ ){
        *t_it = (*oi_it)->transformation_; 
        *sil_it = (*oi_it)->state_ids_; 
      }
      
      op_state_name_canonical_ = "";
      for ( int op_pos : *op_order_ )
        op_state_name_canonical_ += (*op_info_vec_)[op_pos]->op_state_name_; 
     
    } 
    ~MultiOp_Info(){}; 

    std::shared_ptr<std::vector<std::shared_ptr<Op_Info>>> op_info_vec() { return op_info_vec_; }

    //TODO sort out this hack
    char transformation() { assert( op_info_vec_->size() == 1 ); return transformations_->front(); } 

    int op_order( int ii ) { return  (*op_order_)[ii]; }
    std::string op_state_name_canonical() { return op_state_name_canonical_; } 
};

#endif
