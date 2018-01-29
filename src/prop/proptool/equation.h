#ifndef __SRC_PROPTOOL_EQUATION
#define __SRC_PROPTOOL_EQUATION

#include <src/global.h>
#include <src/smith/wicktool/term.h>

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


     virtual void initialize_expression() = 0;
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

     Equation_Init_Value( std::string name,  std::string type, std::shared_ptr<Expression_Init> master_expression, 
                          std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector<int>> >> range_map,
                          std::shared_ptr<std::vector<std::string>> target_indexes, 
                          std::shared_ptr<std::map< std::string, DataType >> factor_map ) :
                          Equation_Init_Base ( name, type, master_expression, range_map ),
                          target_indexes_(target_indexes), factor_map_(factor_map)  {}; 

    ~Equation_Init_Value(){};

     void initialize_expression(); 
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

   
     Equation_Init_LinearRM( std::string name,  std::string type, std::shared_ptr<Expression_Init> master_expression,
                             std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector<int>> >> range_map,
                             std::string target_variable, std::shared_ptr<std::vector<std::string>> target_indexes, 
                             std::shared_ptr<std::map< std::string, DataType >> factor_map ) :
                             Equation_Init_Base ( name, type, master_expression, range_map ),
                             target_variable_(target_variable), target_indexes_(target_indexes),
                             factor_map_(factor_map)  {}; 

    ~Equation_Init_LinearRM(){};

     void initialize_expression(); 
     void build() { std::cout << "Not connected to equation yet" << std::endl;} ; 

}; 
#endif
