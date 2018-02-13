#include <src/prop/proptool/task_translator/equation_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_computers( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_computers " << endl;

  gamma_computer_ = gamma_computer;

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,  
                                                 std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_maps() " << endl;
 
  gamma_data_map_ = gamma_data_map;
  tensop_data_map_ = tensop_data_map;

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::build_expression_computer(){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::build_expression_computer " << endl;

  cout << "equation_->expression_map()->size() = "<< equation_->expression_map()->size() << endl; 
  expression_computer_ =
  make_shared<Expression_Computer::Expression_Computer<DataType>>( gamma_computer_, equation_->expression_map(), range_conversion_map_, tensop_data_map_ );

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
DataType
Equation_Computer_Base<DataType>::get_scalar_result( string result_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 std::cout << "Equation_Computer_Base<DataType>::get_scalar_result" << std::endl;
 DataType x = 1.0 ;
 return x;

 } 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
DataType
Equation_Computer_Base<DataType>::get_scalar_result( string result_name, vector<pair<string, int>>& fixed_idxs,
                                                                         vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 std::cout << "Equation_Computer_Base<DataType>::get_scalar_result" << std::endl;
 DataType x = 1.0 ;
 return x;

 } 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::Tensor_<DataType>> 
Equation_Computer_Base<DataType>::get_tensop( string tensop_name, vector<pair<string, int>>& fixed_idxs,
                                                                  vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop not done" << std::endl;
  std::shared_ptr<SMITH::Tensor_<DataType>> tens = tensop_data_map_->at(tensop_name);

  return tens;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::Tensor_<DataType>> 
Equation_Computer_Base<DataType>::get_tensop( string tensop_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop not done" << std::endl;

  std::shared_ptr<SMITH::Tensor_<DataType>> tens = tensop_data_map_->at(tensop_name);

  return tens;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::MultiTensor_<DataType>>
Equation_Computer_Base<DataType>::get_tensop_vector( string tensop_name, vector<pair<string, int>>& fixed_idxs,
                                                                         vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop not done" << std::endl;

  vector<shared_ptr<Tensor_<DataType>>> tens_list(5); 
  // TODO replace with call to range map in equation for varios input idxs
  vector<int> range = { 1,2,3,4,5}; 
  for ( int  ii : range )  
    tens_list[ii] = tensop_data_map_->at(tensop_name);

  vector<DataType> factor_list( tens_list.size(), 1.0) ; // all ones for now

  shared_ptr<MultiTensor_<DataType>> multitens = make_shared<MultiTensor_<DataType>>( factor_list, tens_list );

  return multitens;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::MultiTensor_<DataType>>
Equation_Computer_Base<DataType>::get_tensop_vector( string tensop_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop not done" << std::endl;

  vector<shared_ptr<Tensor_<DataType>>> tens_list(5); 
  // TODO replace with call to range map in equation for varios input idxs
  vector<int> range = { 1,2,3,4,5}; 
  for ( int  ii : range )  
    tens_list[ii] = tensop_data_map_->at(tensop_name);

  vector<DataType> factor_list( tens_list.size(), 1.0) ; // all ones for now

  shared_ptr<MultiTensor_<DataType>> multitens = make_shared<MultiTensor_<DataType>>( factor_list, tens_list );

  return multitens;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void
Equation_Computer_Base<DataType>::evaluate_term( string expression_name, vector<pair<string, int>>& fixed_idxs,
                                                                               vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::evaluate expression not done" << std::endl;

  expression_computer_->evaluate_expression( expression_name);

  return; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void
Equation_Computer_Base<DataType>::evaluate_term( string expression_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::evaluate expression not done" << std::endl;

  expression_computer_->evaluate_expression( expression_name);

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Value<DataType>::solve_equation(){ /*not really solving anything here*/ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Base<DataType>::Equation_Computer_Value<DataType>::solve_equation() " << endl;

  for ( auto& expr_map_elem : *equation_->expression_map() )
    expression_computer_->evaluate_expression( expr_map_elem.first ) ; 

  return; 
}
////////////////////////////////////////////////////////////////////
template class Equation_Computer_Base<double>;
template class Equation_Computer_Value<double>;
////////////////////////////////////////////////////////////////////
