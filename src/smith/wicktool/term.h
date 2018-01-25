#ifndef __SRC_SMITH_WICKTOOL_TERM_INIT_H
#define __SRC_SMITH_WICKTOOL_TERM_INIT_H

#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/braket.h> // TODO don't like this include... fix equation so it is not necessary

class Op_Init { 
  
  private :
    std::string name_;
    const std::vector<std::string> idxs_;
    const std::vector<int*> idx_ptrs_;
  
  public :
    std::string alg_name_;
  
    Op_Init(std::string base_name, std::vector<std::string>& idxs, std::vector<int*>& idx_ptrs ) :
            name_(base_name), idxs_(idxs), idx_ptrs_(idx_ptrs), alg_name_(base_name) {
              alg_name_ += "_{";
              for (std::string idx : idxs_ ) 
                alg_name_ += idx;
              alg_name_ += "}";
            };

   ~Op_Init(){};

     

     std::string state_specific_name(){
       std::string ssname = name_ ;
       ssname = + "_{";
       for ( auto id_val : idx_ptrs_ )
         ssname += std::to_string(*id_val);
       ssname +="}" ;
       return ssname;
     }

     std::shared_ptr<std::vector<int>> get_idxs() { 
       std::shared_ptr<std::vector<int>> state_idxs_list = std::make_shared<std::vector<int>>(idxs_.size());
       for ( int ii = 0 ; ii != state_idxs_list->size(); ii++ )
         state_idxs_list->at(ii) = *(idx_ptrs_[ii]);
       return state_idxs_list;
     }

}; 


template <typename DataType> 
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

}; 



template <typename DataType> 
class Term_Init { 

  public :
    std::string name_;
    std::string type_;
    std::shared_ptr<std::vector<BraKet_Init<DataType>>> braket_list_;
    std::shared_ptr<std::vector<std::string>> braket_factors_;
    std::shared_ptr<std::map<std::string, int>> index_val_map_;  
   
    std::string alg_name_;

    Term_Init( std::string name, std::string type,
               std::shared_ptr<std::vector<BraKet_Init<DataType>>> braket_list,
               std::shared_ptr<std::vector<std::string>> braket_factors,
               std::shared_ptr<std::map<std::string, int>> index_val_map) :  
               name_(name), type_(type), braket_list_(braket_list), braket_factors_(braket_factors),
               index_val_map_(index_val_map){
             
               alg_name_ = "";
               for ( int ii =0 ; ii != braket_factors_->size(); ii++ )
                 alg_name_ += "(" + braket_factors_->at(ii) + ")" + braket_list->at(ii).name_ + " + ";
               alg_name_.pop_back(); 
               std::cout << "======================= New Term =======================" << std::endl;
               std::cout << alg_name_ << std::endl << std::endl;
               }; 
                                              
    ~Term_Init(){};
}; 


template<typename DataType>
class Expression_Init {

  public : 
    std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init<DataType>>>>> term_list_;
    std::shared_ptr<std::vector<std::map<std::string,std::pair<bool,std::string>>>> term_range_maps_;
   
    std::string name_; 

    Expression_Init( std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init<DataType>>>>> term_list,
                     std::shared_ptr<std::vector<std::map<std::string,std::pair<bool,std::string>>>> term_range_maps ): 
                     term_list_(term_list), term_range_maps_(term_range_maps) {
                       for ( std::pair<std::string, std::shared_ptr<Term_Init<DataType>>> term : *term_list_ ) 
                         name_ += "(" +term.first +")."+ term.second->name_+ "+"; 
                       name_.pop_back();
                       std::cout << "Expression name " << std::endl;
                       std::cout << name_ << std::endl << std::endl;
                     };
    ~Expression_Init(){};
 
};


template<typename DataType>
class Equation_Init {

   public :
     std::string name_;
     std::string type_;
     std::string target_variable_;
     std::map<std::string, std::shared_ptr<std::vector<std::string>>> target_indexes_; // Need a different expression for each one of these.
     std::map<std::string, std::shared_ptr<std::vector<int>>> id_range_map_;   // Need a different expression for each onf of these.
     std::shared_ptr<Expression_Init<DataType>> master_expression_;


     std::shared_ptr<std::vector<std::pair<std::string,std::shared_ptr<Expression_Init<DataType>>>>> expression_list_;
   
     Equation_Init( std::string name,  std::string type, std::string target_variable,
                    std::map<std::string, std::shared_ptr<std::vector<std::string>>> target_indices,
                    std::shared_ptr<Expression_Init<DataType>> master_expression ) :
                    name_(name), type_(type),
                    target_variable_(target_variable), target_indexes_(target_indices), master_expression_(master_expression) {};
     ~Equation_Init(){};


     }; 
 
#endif
