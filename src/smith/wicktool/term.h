#ifndef __SRC_SMITH_WICKTOOL_TERM_INIT_H
#define __SRC_SMITH_WICKTOOL_TERM_INIT_H

#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/braket.h>


template <typename DataType> 
class Term_Init { 

  public :
    std::shared_ptr<std::vector<std::vector< std::pair<std::string, std::vector<std::string>> >>> op_defs_;
    std::vector<std::string> bra_indexes_;
    std::vector<std::string> ket_indexes_;
    std::vector<DataType> factors_;
    std::vector<std::string> types_;
    std::vector<int> all_states_;
  
    std::string name_;

    // Horrific constructor is just building up the name of the term
    Term_Init( std::shared_ptr<std::vector<std::vector< std::pair<std::string, std::vector<std::string>>> >> op_defs,
               std::vector<std::string>& bra_indexes, std::vector<std::string>& ket_indexes,
               std::vector<DataType>& factors, std::vector<std::string>& types) :
               op_defs_(op_defs), bra_indexes_(bra_indexes), ket_indexes_(ket_indexes), factors_(factors), types_(types) { 
        
                 for ( int ii = 0; ii != bra_indexes_.size() ; ii++ ){
                   name_ += +"("+ std::to_string(factors_[ii])+")< " +bra_indexes_[ii] + "| ";
                   for ( std::pair<std::string,std::vector<std::string>> opspec : op_defs_->at(ii) ) {
                     name_ += opspec.first;
                     if (opspec.second.size() != 0 ) {
                        name_ += "_{";  
                        for ( std::string state_id : opspec.second ) 
                          name_ += state_id;
                        name_ += "} ";
                     }
                   }    
		   name_ += "|" + ket_indexes_[ii] + "> + ";
                 }
               }
                


    std::vector<int> merge_states( std::vector<std::vector<int>>& all_st_lists )  { //gross, but it's only ever short
                        std::vector<int> st_list = all_st_lists[0];
                        for ( int ii = 1 ; ii!= all_st_lists.size(); ii++ ){
                          std::vector<int> tmp(st_list.size() + all_st_lists[ii].size());
                          std::copy(st_list.begin(), st_list.end(), tmp.begin() );
                          std::copy(all_st_lists[ii].begin(), all_st_lists[ii].end(), tmp.begin()+st_list.size() );
                          std::sort(tmp.begin(), tmp.end()); 
                          std::vector<int>::iterator it = std::unique(tmp.begin(), tmp.end()); 
                          st_list = std::vector<int>( tmp.begin(), it );
                        }
                        return st_list;
                      }
}; 


template<typename DataType>
class Expression_Init {

  public : 
    std::vector<std::string> BraKet_list_;
    std::shared_ptr<std::vector<std::shared_ptr<Term_Init<DataType>>>> term_list_;
    std::string name_;
   
    Expression_Init( std::shared_ptr<std::vector<std::shared_ptr<Term_Init<DataType>>>> term_list ): 
                     term_list_(term_list){
                       for ( std::shared_ptr<Term_Init<DataType>> term : *term_list_ ) 
                         name_ += term->name_; 
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
    std::map<std::string, std::shared_ptr<std::vector<std::string>>> target_indices_; // Need a different expression for each one of these.
    std::map<std::string, std::shared_ptr<std::vector<std::string>>> summed_indices_; // Need to loop over and add all of these.
    std::map<std::string, std::shared_ptr<std::vector<std::string>>> free_indices_;   // Need a different expression for each onf of these.
    std::shared_ptr<std::vector<std::pair<std::string,std::string>>> expression_list_;

    Equation_Init( std::string name,  std::string type, std::string target_variable,
                   std::map<std::string, std::shared_ptr<std::vector<std::string>>> target_indices,
                   std::map<std::string, std::shared_ptr<std::vector<std::string>>> summed_indices,
                   std::map<std::string, std::shared_ptr<std::vector<std::string>>> free_indices,
                   std::shared_ptr< std::vector<std::pair<std::string,std::string>>> expression_list ) :
                   name_(name), type_(type),
                   target_variable_(target_variable), target_indices_(target_indices), summed_indices_(summed_indices),
                   expression_list_(expression_list) {};
    ~Equation_Init(){};

     std::shared_ptr<std::vector<std::shared_ptr<BraKet<DataType>>>>
     construct_braket_list( std::map<std::string, std::shared_ptr<Expression_Init<DataType>>>& expression_init_map, 
                            std::map<std::string, DataType>& factor_map, std::map<std::string, std::vector<std::string> >& range_map  );


//        for ( std::string exp_name : *expression_list_ ){ 
//          std::shared_ptr<Expression_Init<DataType>> E_init = expression_init_map.at(exp_name); 
//          for ( std::shared_ptr<Term_Init<DataType>> term : *term_list ) {
//            std::vector<std::string> bra_states  = range_map.at(
//
//        op_list_
//
//     } 
                                                     
};
#endif
