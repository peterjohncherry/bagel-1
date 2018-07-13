#include <bagel_config.h>
#include <src/prop/proptool/integrals/moint.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/integrals/moint_computer.h>
#include <src/prop/proptool/debugging_utils.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace bagel; 
using namespace bagel::Tensor_Arithmetic; 
#define __DEBUG_PROPTOOL_MOINT_COMPUTER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void  MOInt_Computer<DataType>::calculate_v2(const vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_v2 IndexRange_ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // TODO for annoying const issue; change arg.
  vector<SMITH::IndexRange>  blocks_buff = blocks;
  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks_buff );

#ifndef __DEBUG_PROPTOOL_MOINT_COMPUTER
  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks_buff );
  vector<int> alt_to_norm_order = { 3, 1, 2, 0 } /* { 2, 0, 3, 1 } */ /* { 2, 1, 3, 0 } */;
  v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order );
  assert(v2_->norm() != 0.0 ) ;
#endif

#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens =   MOInt_Computer<DataType>::calculate_v2_smith() ;
  vector<int> alt_to_norm_order = { 3, 1, 2, 0 };
  v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2_tens, alt_to_norm_order );
  v2_->scale(-1.0);
  assert(v2_->norm() != 0.0 ) ;
#endif

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MOInt_Computer<DataType>::calculate_v2(const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_v2 string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
    blocks[ii] =  *(range_conversion_map_->at(blocks_str[ii])); 
  
  calculate_v2( blocks );

  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is the core ham minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MOInt_Computer<DataType>::calculate_h1( const vector<SMITH::IndexRange>& blocks, bool set_coeffs, bool set_fock ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_h1" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );
  h1_ = one_el_ints.h1();
  
  if ( set_fock ) 
    f1_ = one_el_ints.fock();  

  if ( set_coeffs ) 
    coeffs_ = one_el_ints.coeff();  
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MOInt_Computer<DataType>::calculate_h1( const vector<string>& blocks_str, bool set_coeffs, bool set_fock ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_h1 string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 

  calculate_h1( blocks, set_coeffs, set_fock);

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is core_fock minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MOInt_Computer<DataType>::calculate_fock( const vector<SMITH::IndexRange>& blocks, bool set_coeffs, bool set_h1 ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_fock   : set_coeffs = " << set_coeffs << endl ;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );
  f1_ = one_el_ints.fock();

  if ( set_coeffs ) {  
    coeffs_ = one_el_ints.coeff();  
    got_fock_coeffs_ = true;
  }

  if ( set_h1 ) 
    h1_ = one_el_ints.h1();

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MOInt_Computer<DataType>::calculate_fock( const vector<string>& blocks_str, bool set_coeffs, bool set_h1 ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_fock string ver : set_coeffs = "<< set_coeffs << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 
 
  calculate_fock( blocks, set_coeffs, set_h1); 

  return;
}
/////////////////////////////////////////////////TESTING////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>>  MOInt_Computer<DataType>::calculate_v2_smith() { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_v2_smith IndexRange_ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  SMITH::IndexRange active = *(range_conversion_map_->at("a")); 
  SMITH::IndexRange core = *(range_conversion_map_->at("c")); 
  SMITH::IndexRange virt = *(range_conversion_map_->at("v")); 

  SMITH::IndexRange not_active = core; not_active.merge(virt);
  SMITH::IndexRange not_core = active; not_core.merge(virt); 
  SMITH::IndexRange not_virt = core;   not_virt.merge(active);
   
  SMITH::IndexRange free = core; free.merge(active); free.merge(virt);
   
  // TODO for annoying const issue; change arg.
  vector<SMITH::IndexRange>  blocks = { free, free, free, free };
  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );

  auto v2_tens = v2.tensor();
  
  //////////////////////////////////////////////////////////////////////
  vector<SMITH::IndexRange> keep_ranges = { virt, virt, core, core };
  //////////////////////////////////////////////////////////////////////

  {
    vector<SMITH::IndexRange>  keep_ranges_smith_order = { keep_ranges[3], keep_ranges[1], keep_ranges[2], keep_ranges[0] };
    auto v2_keep = Tensor_Arithmetic_Utils::get_sub_tensor( v2_tens, keep_ranges_smith_order  );
    v2_tens->zero();
    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( v2_keep, v2_tens ); 
  }
  return v2_tens;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_test_tensor( const vector<SMITH::IndexRange>& index_ranges ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_test_tensor string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  return  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_tensor_row_major( index_ranges );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>                                     
shared_ptr<SMITH::Tensor_<DataType>>
MOInt_Computer<DataType>::build_s_test_tensor( const vector<int>& ordering) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::build_s_test_tensor" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  SMITH::IndexRange active = *(range_conversion_map_->at("a")); 
  SMITH::IndexRange core = *(range_conversion_map_->at("c")); 
  SMITH::IndexRange virt = *(range_conversion_map_->at("v")); 


  SMITH::IndexRange not_active = core; not_active.merge(virt);
  SMITH::IndexRange not_core = active; not_core.merge(virt); 
  SMITH::IndexRange not_virt = core; not_virt.merge(active);

  // Smith order : {3, 1, 2, 0}  
  vector<SMITH::IndexRange> full_block = { not_core, not_core, not_virt, not_virt };
  // vector<SMITH::IndexRange> full_block = { not_virt, not_virt, not_core, not_core };
  shared_ptr<SMITH::Tensor_<DataType>> s_test_tensor =  make_shared<SMITH::Tensor_<DataType>>( WickUtils::reorder_vector( ordering, full_block ) ); 
  s_test_tensor->allocate();
  s_test_tensor->zero();

//  {
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, core, active }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(1.0) );
//  }{
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, active, core }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(-1.0) );
//  }{
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, active, active }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(0.0) );
//  }
  {
  vector<SMITH::IndexRange> non_zero_block = { virt, virt, core, core }; 
  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(1.0) );
  }
  cout << "  s_test_tensor->norm() = " ; cout.flush(); cout <<  s_test_tensor->norm() <<  endl;
 
  return s_test_tensor; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_test_tensor( const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_test_tensor string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 

  return  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_tensor_row_major( blocks );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class MOInt_Computer<double>;
template class MOInt_Computer<std::complex<double>>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
