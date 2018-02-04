#include <bagel_config.h>
#include  <src/prop/proptool/task_translator/system_computer.h>


using namespace std;
using namespace bagel;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
System_Computer::System_Computer<DataType>::System_Computer(
                    shared_ptr<System_Info<DataType>> system_info,
                    shared_ptr<MOInt_Computer<double>> moint_computer, 
                    shared_ptr<map< string, shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                    shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer ):
                    system_info_(system_info), moint_computer_(moint_computer), range_conversion_map_(range_conversion_map),
                    b_gamma_computer_(b_gamma_computer) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Computer::System_Computer<DataType>::System_Computer()" << endl; 
  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();
  gamma_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<double>>>>();

  b_gamma_computer_->set_maps( range_conversion_map_, system_info_->Gamma_map, gamma_data_map_,
                               sigma_data_map_, civec_data_map_ );

  calculate_mo_integrals();

 cout << "should either set or initialize maps here " << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void System_Computer::System_Computer<DataType>::build_equation_computer(std::string equation_name ){ 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " void System_Computer::System_Computer<DataType>::build_equation_computer" << endl; 

  shared_ptr<Equation_Base<DataType>> equation_basic = system_info_->equation_map_->at(equation_name); 

  if ( equation_basic->type()  ==  "Value" ) { 
    //shared_ptr<Equation_Value<DataType>> equation_val = dynamic_pointer_cast<Equation_Value<DataType>>(equation_basic); 
    shared_ptr<Equation_Computer_Value<DataType>> equation_computer = make_shared<Equation_Computer_Value<DataType>>( equation_basic, range_conversion_map_ ); 
    equation_computer->set_computers( b_gamma_computer_ ) ;
    equation_computer->set_maps( gamma_data_map_, tensop_data_map_ );
    equation_computer->build_expression_computer();
    equation_computer->solve_equation(); 
  }
//  dynamic_pointer_cast<Equation_Base<DataType>>(equation(make_shared<Equation_Computer_Value<DataType>>( equation, range_conversion_map_ ))); 
//    equation_computer = dynamic_pointer_cast<Equation_Computer_Base<DataType>>(make_shared<Equation_Computer_Value<DataType>>( equation, range_conversion_map_ )); 
   
//  } else if ( equation->type()  ==  "LinearRM" ) { 
 //   equation_computer = dynamic_pointer_cast<Equation_Computer_Base<DataType>>(make_shared<Equation_Computer_LinearRM<DataType>>( equation, range_conversion_map_ )); 

//  } else {  
//    throw std::logic_error( "this type of equation has not been implemented yet! Aborting!!" );  
//  }

  return;
} 

////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void System_Computer::System_Computer<DataType>::build_expression_computer( std::string expression_name ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " void System_Computer::System_Computer<DataType>::build_expression_computer" << endl; 

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void System_Computer::System_Computer<DataType>::build_tensop( std::string tensop_name ){ 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " void System_Computer::System_Computer<DataType>::build_tensop" << endl; 

  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void System_Computer::System_Computer<DataType>::calculate_mo_integrals() {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Computer::System_Computer::calculate_mo_integrals()" << endl;

  cout << "getting mo integrals " <<  endl; 
  vector<string> test_ranges4 = { "notcor", "notcor", "notvir", "notvir" }; 
  vector<string> test_ranges2 = { "free", "free" }; 
  h1_  =  moint_computer_->get_h1( test_ranges2, true ) ;
  v2_  =  moint_computer_->get_v2( test_ranges4 ) ;
  cout << " new_coeffs  v2->norm() = " << v2_->norm() << endl; 

  tensop_data_map_->emplace( "H" , v2_ );  
  tensop_data_map_->emplace( "f" , h1_ );  
   
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Computer::System_Computer<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
