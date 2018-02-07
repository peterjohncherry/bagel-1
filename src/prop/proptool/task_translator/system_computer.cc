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

  //TODO Skipping for now
  // TODO Should only set the tensors on an expression by expression basis
  // Currently gets all necessary tensor blocks at once...
  //  for ( auto& expr_map_loc : *equation_basic->expression_map() )
  //  get_necessary_tensor_blocks(expr_map_loc.second);

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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void System_Computer::System_Computer<DataType>::get_necessary_tensor_blocks( shared_ptr<Expression<DataType>> expression ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Expression_Computer::Expression_Computer<DataType>::set_gamma_maps" << endl;

  // Loop through A-tensors needed for this gamma
  for ( auto& gamma_contribs : *expression->G_to_A_map_ ){
    for ( auto& A_contrib_list : *gamma_contribs.second ){
      auto& A_contrib  = A_contrib_list.second; 
      for ( auto& CTP : *( expression->CMTP_map()->at(A_contrib.name())->get_CTP_vec() ) )
        if ( tensop_data_map_->find(CTP->myname() ) == tensop_data_map_->end() )
          get_tensor_block( CTP->myname(), CTP->full_id_ranges() );
    }
  }
  return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void System_Computer::System_Computer<DataType>::get_tensor_block( string tensor_block_name, shared_ptr<vector<string>> idx_ranges ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "System_Computer::System_Computer::get_tensor_block" << endl;

  shared_ptr<vector<SMITH::IndexRange>> smith_idx_ranges = make_shared<vector<SMITH::IndexRange>>(idx_ranges->size());
  vector<SMITH::IndexRange>::iterator sir_it = smith_idx_ranges->begin();
  for ( string& range_name :  *idx_ranges )
    *sir_it++ = *(range_conversion_map_->at(range_name));
 
  if (tensor_block_name[0] == 'H' ) {
  //  tensop_data_map_->emplace( tensor_block_name, moint_computer_->get_v2( *idx_ranges ) );

  shared_ptr<SMITH::Tensor_<double>> H_block = Tensor_Arithmetic_Utils::get_sub_tensor( tensop_data_map_->at("H") , *smith_idx_ranges ) ;

  } else {
   
    DataType one = (DataType)(1.0); //TODO find a better way;
    tensop_data_map_->emplace( "T" , Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_Tensor( smith_idx_ranges, one ) );
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
  vector<string> act4 = { "act", "act", "act", "act" };
  vector<string> free2 = { "free", "free" };

  h1_  =  moint_computer_->get_h1( free2, true ) ;
  v2_  =  moint_computer_->get_v2( free4 ) ;
  cout << " new_coeffs  v2->norm() = " << v2_->norm() << endl;

  tensop_data_map_->emplace( "H" , v2_ );
  tensop_data_map_->emplace( "f" , h1_ );

//  tensop_data_map_->emplace( "T" , moint_computer_->get_test_tensor( free4 ) );
  //tensop_data_map_->emplace( "X" , moint_computer_->get_test_tensor( free4 ) );
  DataType one = (DataType)(1.0); //TODO find a better way;
  SMITH::IndexRange fs = *(range_conversion_map_->at("free"));
  SMITH::IndexRange nvs = *(range_conversion_map_->at("cor"));
  nvs.merge(*(range_conversion_map_->at("act")));

  SMITH::IndexRange ncs = *(range_conversion_map_->at("act"));
  ncs.merge(*(range_conversion_map_->at("virt")));

  shared_ptr<vector<SMITH::IndexRange>> fs4 = make_shared<vector<SMITH::IndexRange>>(vector<SMITH::IndexRange> { ncs, ncs, nvs, nvs } );   
  tensop_data_map_->emplace( "X" , Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_Tensor( fs4, one ) );

  cout <<"X->norm() = " << tensop_data_map_->at("X")->norm() << endl; 

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class System_Computer::System_Computer<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
