#include <bagel_config.h>
#include  <src/prop/proptool/task_translator/system_computer.h>


using namespace std;
using namespace bagel;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
System_Computer::System_Computer<DataType>::System_Computer(shared_ptr<System_Info<DataType>> system_info,
                                                            shared_ptr<MOInt_Computer<DataType>> moint_computer,
                                                            shared_ptr<map< string, shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                                                            shared_ptr<B_Gamma_Computer::B_Gamma_Computer<DataType>> b_gamma_computer ):
                                                            system_info_(system_info), moint_computer_(moint_computer),
                                                            range_conversion_map_(range_conversion_map),
                                                            b_gamma_computer_(b_gamma_computer) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL
cout << "System_Computer::System_Computer<DataType>::System_Computer()" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

  gamma_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

  b_gamma_computer_->set_maps( range_conversion_map_, system_info_->Gamma_map, gamma_data_map_,
                               sigma_data_map_, civec_data_map_ );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void System_Computer::System_Computer<DataType>::build_equation_computer(std::string equation_name ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_SYSTEM_COMPUTER
cout << "System_Computer<DataType>::build_equation_computer : " << equation_name << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Equation_Base<DataType>> equation_basic = system_info_->equation_map_->at(equation_name);

  if ( equation_basic->type()  ==  "Value" ) {
    shared_ptr<Equation_Computer_Value<DataType>> equation_computer = make_shared<Equation_Computer_Value<DataType>>( equation_basic, range_conversion_map_ );
    equation_computer->set_computers( b_gamma_computer_, moint_computer_ ) ;
    equation_computer->set_maps( gamma_data_map_, tensop_data_map_ );
    equation_computer->build_expression_computer();
    equation_computer->solve_equation();
  
  } else if ( equation_basic->type()  ==  "LinearRM" ) {
    auto equation_linearrm = dynamic_pointer_cast<Equation_LinearRM<DataType>>( equation_basic ) ;
    if (!equation_linearrm) { throw logic_error("casting of sptr<Equation_base> to sptr<Equation_LinearRM>  failed") ; } 

    shared_ptr<Equation_Computer_LinearRM<DataType>> equation_computer = make_shared<Equation_Computer_LinearRM<DataType>>( equation_linearrm, range_conversion_map_ );
    equation_computer->set_computers( b_gamma_computer_, moint_computer_ ) ;
    equation_computer->set_maps( gamma_data_map_, tensop_data_map_ );
    equation_computer->build_expression_computer(); // TODO WRONG NAME AND WRONG FUNCTIONALITY
    equation_computer->solve_equation();
   
  } else {
    throw std::logic_error( "this type of equation has not been implemented yet! Aborting!!" );
  }
 
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Computer::System_Computer<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
