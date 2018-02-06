#ifndef __SRC_PROPTOOL_EQUATION
#define __SRC_PROPTOOL_EQUATION

#include <src/global.h>
#include <src/prop/proptool/initialization/op_bk_term_expr_init.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
// Equation_Init constructs the expressions necessary for evaluation of the equation
// specified in the user input.
// The "master expression" is the expression with variables for indexes.
// This is used to generate a set of expressions for which the numbers are specified.
// These expressions are defined in by a list of terms, which are themselves a list of
// brakets. 
// The State_numbers in the BraKets are not variables, but numbers, which are determined
// by looping through the vectors in the range map (the loop structure is determined from
// the summations specified in the input). 
// The Terms and BraKets are defined seperately and kept, as the same BraKets may 
// appear in different expressions, and so can be reused.
class Equation_Init_Base {

   public :

     std::string name_;
     std::string type_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<int>>>> range_map_; 
     std::shared_ptr<Expression_Init> master_expression_;

     Equation_Init_Base( std::string name, std::string type, 
                         std::shared_ptr<Expression_Init> master_expression,
                         std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<int>>>> range_map ) :
                         name_(name), type_(type),  master_expression_(master_expression),
                         range_map_(range_map) {};

     ~Equation_Init_Base(){};


     virtual void initialize_expressions() = 0;
     virtual void build() = 0;

}; 
 
// Generates an Equation object to evaluate all f_ij 
// f is the master expression
// i and j range over all values specified by target indexes
template<typename DataType>
class Equation_Init_Value : public Equation_Init_Base {

   public :
   
     DataType factor_;
     std::shared_ptr<std::vector<std::string>> target_indexes_;                // Need a different expression for each one of these.
     std::shared_ptr<std::map< std::string, DataType >> factor_map_; 
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_;
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_;

     Equation_Init_Value( std::string name,  std::string type, std::shared_ptr<Expression_Init> master_expression, 
                          std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector<int>> >> range_map,
                          std::shared_ptr<std::vector<std::string>> target_indexes, 
                          std::shared_ptr<std::map< std::string, DataType >> factor_map ) :
                          Equation_Init_Base ( name, type, master_expression, range_map ),
                          target_indexes_(target_indexes), factor_map_(factor_map)  {

                          expression_term_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>>();
                          term_braket_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>();

                          }; 

    ~Equation_Init_Value(){};

     void initialize_expressions(); 
     void build() { std::cout << "Not connected to equation yet" << std::endl;} ; 


}; 

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
     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_;

   
     Equation_Init_LinearRM( std::string name,  std::string type, std::shared_ptr<Expression_Init> master_expression,
                             std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector<int>> >> range_map,
                             std::string target_variable, std::shared_ptr<std::vector<std::string>> target_indexes, 
                             std::shared_ptr<std::map< std::string, DataType >> factor_map ) :
                             Equation_Init_Base ( name, type, master_expression, range_map ),
                             target_variable_(target_variable), target_indexes_(target_indexes),
                             factor_map_(factor_map)  {

                             expression_term_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>>();
                             term_braket_map_ = std::make_shared<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>();

}; 

    ~Equation_Init_LinearRM(){};

     void initialize_expressions(); 
     void build() { std::cout << "Not connected to equation yet" << std::endl;} ; 

}; 
#endif
