#ifndef __SRC_SMITH_WICKTOOL_TERM_INIT_H
#define __SRC_SMITH_WICKTOOL_TERM_INIT_H

#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/braket.h>

class Op_Init { 
  
  private :
    const std::vector<std::string> idxs_;
    const std::shared_ptr<std::vector<int*>> idx_ptrs_;
  
  public :
    std::string name_;
    std::string alg_name_;
    int state_dep_ = 0;
     
    Op_Init(std::string base_name, std::vector<std::string>& idxs, std::shared_ptr<std::vector<int*>> idx_ptrs ) :
            name_(base_name), idxs_(idxs), idx_ptrs_(idx_ptrs), alg_name_(base_name) {
              if (idxs.front() != "none" ) { 
                alg_name_ += "_{";
                for (std::string idx : idxs_ ) 
                  alg_name_ += idx;
                alg_name_ += "}";
                state_dep_ = idxs_.size(); 
              } else { 
                state_dep_ = 0; 
              }
            };

   ~Op_Init(){};

     int num_idxs(){ return idxs_.size(); }   

     std::string state_specific_name(){
       std::string ssname = name_ ;
       ssname = + "_{";
       for ( auto id_val : *idx_ptrs_ )
         ssname += std::to_string(*id_val);
       ssname +="}" ;
       return ssname;
     }

     std::shared_ptr<std::vector<int>> get_idxs() { 
       std::shared_ptr<std::vector<int>> state_idxs_list = std::make_shared<std::vector<int>>(idxs_.size());
       if (idxs_.front() != "none" ) {
         for ( int ii = 0 ; ii != state_idxs_list->size(); ii++ )
         state_idxs_list->at(ii) = *(idx_ptrs_->at(ii));
       }  
       return state_idxs_list;
     }

     void get_op_idxs(std::vector<int>& state_idxs_list ) { // to avoid incredible number of shared_ptrs 
       if (idxs_.front() != "none" ) {
         for ( int ii = 0 ; ii != state_idxs_list.size(); ii++ ) 
           state_idxs_list[ii] = *(idx_ptrs_->at(ii));
       }  
       return;
     }

}; 


class BraKet_Init { 

  public :

   std::shared_ptr<std::vector<Op_Init>> op_list_;

   std::string bra_index_;
   int* bra_index_ptr_;

   std::string ket_index_;
   int* ket_index_ptr_;

   std::string name_;

   BraKet_Init( std::shared_ptr<std::vector<Op_Init>> op_list, 
                std::string bra_index, int* bra_index_ptr, std::string ket_index, int* ket_index_ptr ) :
                op_list_(op_list), bra_index_(bra_index), bra_index_ptr_(bra_index_ptr),
                ket_index_(ket_index),  ket_index_ptr_(ket_index_ptr) {
               
                  name_ = "<" + bra_index + "|";
                  for ( Op_Init op : *op_list )
                    name_ += op.alg_name_;    
                  
                  name_+= "|" + ket_index_+ ">" ;
               };

   ~BraKet_Init(){};

   std::shared_ptr<std::vector<std::vector<int>>> get_op_idxs_list() { 
     std::shared_ptr<std::vector<std::vector<int>>> op_idxs_list = std::make_shared<std::vector<std::vector<int>>>(op_list_->size());
     for ( int ii = 0 ; ii != op_idxs_list->size(); ii++ ){
       op_idxs_list->at(ii) = std::vector<int>( op_list_->at(ii).num_idxs() );
       op_list_->at(ii).get_op_idxs( op_idxs_list->at(ii) );
     } 
     return op_idxs_list;
   }

   void get_op_idxs_list( std::vector<std::vector<int>>& op_idxs_list ) { 
     for ( int ii = 0 ; ii != op_idxs_list.size(); ii++ ){
       op_idxs_list[ii] = std::vector<int>( op_list_->at(ii).num_idxs() );
       op_list_->at(ii).get_op_idxs( op_idxs_list[ii] );
     } 
     return;
   }

   int bra_index(){ return *bra_index_ptr_; }
   int ket_index(){ return *ket_index_ptr_; }
   

}; 

class Term_Init { 

  public :
    std::string name_;
    std::string type_;
    std::shared_ptr<std::vector<BraKet_Init>> braket_list_;
    std::shared_ptr<std::vector<std::string>> braket_factors_;
    std::shared_ptr<std::map<std::string, int>> idx_val_map_;  
   
    std::string alg_name_;

    Term_Init( std::string name, std::string type,
               std::shared_ptr<std::vector<BraKet_Init>> braket_list,
               std::shared_ptr<std::vector<std::string>> braket_factors,
               std::shared_ptr<std::map<std::string, int>> idx_val_map) :  
               name_(name), type_(type), braket_list_(braket_list), braket_factors_(braket_factors),
               idx_val_map_(idx_val_map){
             
               alg_name_ = "";
               for ( int ii =0 ; ii != braket_factors_->size(); ii++ )
                 alg_name_ += "(" + braket_factors_->at(ii) + ")" + braket_list->at(ii).name_ + " + ";
               alg_name_.pop_back(); 
               std::cout << "======================= New Term =======================" << std::endl;
               std::cout << alg_name_ << std::endl << std::endl;
               }; 

                                               

    ~Term_Init(){};
}; 

class Expression_Init {

  public : 
    std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init>>>> term_list_;
    std::shared_ptr<std::vector<std::shared_ptr<std::map<std::string,std::pair<bool,std::string>>>>> term_range_maps_;
   
    std::string name_; 

    Expression_Init( std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init>>>> term_list,
                     std::shared_ptr<std::vector<std::shared_ptr<std::map<std::string,std::pair<bool,std::string>>>>> term_range_maps ): 
                     term_list_(term_list), term_range_maps_(term_range_maps) {
                       for ( std::pair<std::string, std::shared_ptr<Term_Init>> term : *term_list_ ) 
                         name_ += "(" +term.first +")."+ term.second->name_+ "+"; 
                       name_.pop_back();
                       std::cout << "Expression name " << std::endl;
                       std::cout << name_ << std::endl << std::endl;
                     };
    ~Expression_Init(){};
 
};
#endif
