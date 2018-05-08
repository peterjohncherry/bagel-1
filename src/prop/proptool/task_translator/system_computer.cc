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
  cout << "System_Computer::System_Computer<DataType>::System_Computer()" << endl;

  sigma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();
  civec_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();
  gamma_data_map_  = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

  tensop_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

  gamma_data_map_ = make_shared<map< string, shared_ptr<SMITH::Tensor_<DataType>>>>();

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
    shared_ptr<Equation_Computer_Value<DataType>> equation_computer = make_shared<Equation_Computer_Value<DataType>>( equation_basic, range_conversion_map_ );
    equation_computer->set_computers( b_gamma_computer_ ) ;
    equation_computer->set_maps( gamma_data_map_, tensop_data_map_ );
    equation_computer->build_expression_computer();
    equation_computer->solve_equation();
  
  } else if ( equation_basic->type()  ==  "LinearRM" ) {
    cout << " building linear RM computer " << endl;
    auto equation_linearrm = dynamic_pointer_cast<Equation_LinearRM<DataType>>( equation_basic ) ;
    if (!equation_linearrm) { throw logic_error("casting of sptr<Equation_base> to sptr<Equation_LinearRM>  failed") ; } 

    shared_ptr<Equation_Computer_LinearRM<DataType>> equation_computer = make_shared<Equation_Computer_LinearRM<DataType>>( equation_linearrm, range_conversion_map_ );
    equation_computer->set_computers( b_gamma_computer_ ) ;
    equation_computer->set_maps( gamma_data_map_, tensop_data_map_ );
    equation_computer->build_expression_computer(); // TODO WRONG NAME AND WRONG FUNCTIONALITY
    equation_computer->solve_equation();
   
  } else {
  
    throw std::logic_error( "this type of equation has not been implemented yet! Aborting!!" );
  }
 
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
  vector<string> free4 = { "free", "free", "free", "free" };
  vector<string> free2 = { "free", "free" };

  f1_ =  moint_computer_->get_fock( free2, false ) ;
  h1_  =  moint_computer_->get_h1( free2, true ) ;
  v2_  =  moint_computer_->get_v2( free4 ) ;

  //TODO Bad test for now; put the state symmetry into gammagen
  tensop_data_map_->emplace( "H_{00}" , v2_ );
  tensop_data_map_->emplace( "H_{01}" , v2_ );
  tensop_data_map_->emplace( "H_{10}" , v2_ );
  tensop_data_map_->emplace( "H_{11}" , v2_ );
  tensop_data_map_->emplace( "h" , h1_ );
  tensop_data_map_->emplace( "f" , f1_ );

  //  tensop_data_map_->emplace( "T" , moint_computer_->get_test_tensor( free4 ) );
  //tensop_data_map_->emplace( "X" , moint_computer_->get_test_tensor( free4 ) );
  DataType one = (DataType)(1.0); //TODO find a better way;
  SMITH::IndexRange fs = *(range_conversion_map_->at("free"));
  SMITH::IndexRange nvs = *(range_conversion_map_->at("c"));
  nvs.merge(*(range_conversion_map_->at("a")));

  SMITH::IndexRange ncs = *(range_conversion_map_->at("a"));
  ncs.merge(*(range_conversion_map_->at("v")));

  shared_ptr<vector<SMITH::IndexRange>> fs4 = make_shared<vector<SMITH::IndexRange>>(vector<SMITH::IndexRange> { ncs, ncs, nvs, nvs } );   
  shared_ptr<SMITH::Tensor_<DataType>> XTens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_Tensor( fs4, one ); 
  tensop_data_map_->emplace( "X" , XTens );

  cout <<"X->norm() = "; cout.flush() ; cout << tensop_data_map_->at("X")->norm() << endl; 

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Computer::System_Computer<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
