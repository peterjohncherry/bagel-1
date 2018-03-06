#ifndef __SRC_PROPTOOL_ALGMANIP_EQUATION_INIT_LRM
#define __SRC_PROPTOOL_ALGMANIP_EQUATION_INIT_LRM

#include<src/prop/proptool/initialization/equation_init.h>

// Will solve f[T_{ij}] = 0  for T_{ij}
// f is the master_expression.
// T is the target variable. 
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Init_LinearRM : public Equation_Init_Base {

   public :

     std::string target_variable_;
     std::shared_ptr<std::vector<std::string>> target_indexes_;                // Need a different expression for each one of these.
     std::shared_ptr<std::map< std::string, DataType >> factor_map_; 

     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_;

     std::shared_ptr<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>, 
                               std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_state_spec_;

     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_;

     std::shared_ptr<std::map<std::pair< std::string, std::vector<std::pair<std::string, int>>>, 
                                         std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_state_spec_;
   
     Equation_Init_LinearRM( std::string name,  std::string type, std::shared_ptr<std::vector<std::shared_ptr<Expression_Init>>> master_expressions,
                             std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector<int>> >> range_map,
                             std::string target_variable, std::shared_ptr<std::vector<std::string>> target_indexes, 
                             std::shared_ptr<std::map< std::string, DataType >> factor_map ) :
                             Equation_Init_Base ( name, type, master_expressions, range_map ),
                             target_variable_(target_variable), target_indexes_(target_indexes),
                             factor_map_(factor_map)  {

                             term_braket_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>();

                             term_braket_map_state_spec_ = std::make_shared<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>,
                                                                                     std::shared_ptr<std::vector<BraKet<DataType>>>>>();

                             expression_term_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>>();

                             expression_term_map_state_spec_ = std::make_shared<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>, 
                                                                                         std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>>();
                             }; 

    ~Equation_Init_LinearRM(){};

     void initialize_expressions(); 
     void initialize_all_terms(); 
     void build() { std::cout << "Not connected to equation yet" << std::endl;} ; 

}; 
#endif
