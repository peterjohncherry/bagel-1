#ifndef __SRC_SMITH_WICKTOOL_TERM_INIT_H
#define __SRC_SMITH_WICKTOOL_TERM_INIT_H

#include <src/smith/wicktool/wickutils.h>

template <typename DataType> 
  class Term_Init { 
  
    public :
      std::vector<std::vector<std::string>> op_names_;
      std::vector<std::vector<int>> bra_states_;
      std::vector<std::vector<int>> ket_states_;
      std::string bra_index_;
      std::string ket_index_;
      std::vector<DataType> factors_;
      std::vector<std::string> types_;
      std::vector<int> bra_states_merged_;
      std::vector<int> ket_states_merged_;
      std::vector<int> all_states_;
  
      Term_Init( std::vector<std::vector<std::string>>& op_names, std::vector<std::vector<int>>& bra_states, std::vector<std::vector<int>>& ket_states,
                 std::vector<DataType>& factors, std::vector<std::string> types) :
                 op_names_(op_names), bra_states_(bra_states), ket_states_(ket_states), factors_(factors), types_(types) {
  
                  bra_states_merged_ = merge_states(bra_states_);
                  ket_states_merged_ = merge_states(ket_states_);
  
                  std::vector<std::vector<int>> tmp  = { bra_states_merged_, ket_states_merged_ };
                  all_states_ = merge_states(tmp);
  
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
class Equation_Init {

      public : 

        std::string target_variable_;
        std::string target_ranges_;
        std::string summed_ranges_;
        std::string method_;
        std::vector<std::string> expression_list_;
        
        Equation_Init( std::string target_variable, std::string target_ranges,  std::string summed_ranges,  std::string method, std::vector<std::string>& expression_list ) :
                       target_variable_(target_variable), target_ranges_(target_ranges), summed_ranges_(summed_ranges), method_(method), 
                       expression_list_(expression_list) {}; 
        ~Equation_Init() {} ;
 
};
#endif
