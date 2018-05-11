#include <src/prop/proptool/task_translator/equation_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

//TODO Several routines here should be changed so that they use the OpInfo or MultiOpInfo class. 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_computers( std::shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> gamma_computer, 
                                                      std::shared_ptr<MOInt_Computer<DataType>> moint_computer) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_computers " << endl;

  gamma_computer_ = gamma_computer;
  moint_computer_ = moint_computer;

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Computer_Base<DataType>::set_maps( shared_ptr<map< string, shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,  
                                                 shared_ptr<map< string, shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Value<DataType>::set_maps() " << endl;
 
  gamma_data_map_ = gamma_data_map;
  system_tensop_data_map_ = tensop_data_map;
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
  make_shared<Expression_Computer::Expression_Computer<DataType>>( gamma_computer_, equation_->expression_map(), range_conversion_map_, tensop_data_map_,
                                                                   moint_computer_ );

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
DataType
Equation_Computer_Base<DataType>::get_scalar_result( string result_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 std::cout << "Equation_Computer_Base<DataType>::get_scalar_result (name, fixed) "  << result_name << std::endl;
 throw logic_error("not implemented");
 DataType x = 1.0 ;
 return x;

} 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
DataType
Equation_Computer_Base<DataType>::get_scalar_result( string result_name, vector<pair<string, int>>& fixed_idxs,
                                                                         vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 std::cout << "Equation_Computer_Base<DataType>::get_scalar_result (name, fixed, summed)" << result_name <<  std::endl;
 throw logic_error("not implemented");
 DataType x = 1.0 ;
 return x;

 } 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::Tensor_<DataType>> 
Equation_Computer_Base<DataType>::get_tensop( string tensop_name, vector<pair<string, int>>& fixed_idxs,
                                                                  vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop( name, fixed, summed )" << std::endl;
  assert( tensop_data_map_->find(tensop_name) != tensop_data_map_->end() ) ;
  std::shared_ptr<SMITH::Tensor_<DataType>> tens = tensop_data_map_->at(tensop_name);

  return tens;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::Tensor_<DataType>> 
Equation_Computer_Base<DataType>::get_tensop( string tensop_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop (name, fixed)" << std::endl;

  assert( tensop_data_map_->find(tensop_name) != tensop_data_map_->end() ) ;
  std::shared_ptr<SMITH::Tensor_<DataType>> tens = tensop_data_map_->at(tensop_name);

  return tens;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::MultiTensor_<DataType>>
Equation_Computer_Base<DataType>::get_tensop_vector( string tensop_name, vector<pair<string, int>>& fixed_idxs,
                                                                         vector<pair<string, int>>& summed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer_Base<DataType>::get_tensop_vector ( name, fixed, summed) : tensop_name = "; cout.flush(); cout << tensop_name << endl;

  vector<shared_ptr<Tensor_<DataType>>> tens_list(5); 
  // TODO replace with call to range map in equation for varios input idxs
  vector<int> range = { 1,2,3,4,5}; 
  for ( int  ii : range ) { 
    assert( tensop_data_map_->find(tensop_name) != tensop_data_map_->end() ) ;
    tens_list[ii] = tensop_data_map_->at(tensop_name);
  }
  vector<DataType> factor_list( tens_list.size(), 1.0) ; // all ones for now

  shared_ptr<MultiTensor_<DataType>> multitens = make_shared<MultiTensor_<DataType>>( factor_list, tens_list );

  return multitens;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
std::shared_ptr<SMITH::MultiTensor_<DataType>>
Equation_Computer_Base<DataType>::get_tensop_vector( string tensop_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::get_tensop_vector (name, fixed) : tensop_name = "; cout.flush(); cout  << tensop_name << endl;  

  vector<shared_ptr<Tensor_<DataType>>> tens_list(5); 
  // TODO replace with call to range map in equation for varios input idxs
  vector<int> range = { 1,2,3,4,5 }; 
  for ( int  ii : range ) { 
    assert( tensop_data_map_->find(tensop_name) != tensop_data_map_->end() ) ;
    tens_list[ii] = tensop_data_map_->at(tensop_name);
  }

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
  std::cout << "Equation_Computer_Base<DataType>::evaluate_term A : " << expression_name << std::endl;

  expression_computer_->evaluate_expression( expression_name);

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void
Equation_Computer_Base<DataType>::evaluate_term( string expression_name, vector<pair<string, int>>& fixed_idxs ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Equation_Computer_Base<DataType>::evaluate_term B : " << expression_name << std::endl;

  cout << "expression_name = " << expression_name  << " = [ ";  cout.flush();
  for ( auto& xx : fixed_idxs ) {
    cout << "(" << xx.first << ","<< xx.second << ") "; cout.flush();
  } 
  cout << "]" << endl;

  shared_ptr<Expression<DataType>> expr =  equation_->get_term( expression_name, fixed_idxs);

  cout << "got expression : expr->name() = " << expr->name() <<endl;  

  expression_computer_->evaluate_expression(equation_->get_term( expression_name, fixed_idxs));

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
