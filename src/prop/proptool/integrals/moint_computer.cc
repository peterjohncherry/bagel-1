#include <bagel_config.h>
#include <src/prop/proptool/integrals/moint.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/integrals/moint_computer.h>
#include <src/prop/proptool/debugging_utils.h>

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
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens =   MOInt_Computer<DataType>::calculate_v2_smith( blocks_buff ) ;
  vector<int> alt_to_norm_order = { 3, 1, 2, 0 } /* { 2, 0, 3, 1 } */ /*{ 2, 1, 3, 0 } */ ;
  v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2_tens, alt_to_norm_order );
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::calculate_v2_smith(const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_v2_smith string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
    blocks[ii] =  *(range_conversion_map_->at(blocks_str[ii])); 
 
  return calculate_v2_smith( blocks );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>>  MOInt_Computer<DataType>::calculate_v2_smith(const vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::calculate_v2_smith IndexRange_ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // TODO for annoying const issue; change arg.
  vector<SMITH::IndexRange>  blocks_buff = blocks;
  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks_buff );

  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = v2.tensor();
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
  
    shared_ptr<SMITH::Tensor_<DataType>> v2_new = get_test_tensor( blocks_buff );
    SMITH::IndexRange active_ = *(range_conversion_map_->at("a")); 
    SMITH::IndexRange core_ = *(range_conversion_map_->at("c")); 
    SMITH::IndexRange virt_ = *(range_conversion_map_->at("v")); 

    // Block to keep 
    vector<SMITH::IndexRange> keep_ranges = {active_, virt_, active_, virt_};

    shared_ptr<SMITH::Tensor_<DataType>> v2_tens_keep = Tensor_Arithmetic_Utils::get_sub_tensor( v2_new, keep_ranges );

    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( v2_tens_keep, 0.5 );

    cout << endl << endl;
    cout << "======================= v2_->norm() ========================" << endl;
    cout  << " v2_->norm() = " << v2_new->norm() << endl;
    cout  << " v2_tens_keep->norm() = " << v2_tens_keep->norm() << endl;
    v2_new->zero();
    cout << "POST ZEROING" << endl;
    cout << " v2_->norm() = " << v2_tens_keep->norm() <<  endl;
    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( v2_tens_keep, v2_new );
    cout << "POST PUTTING" << endl;
    cout  << " v2.tensor()->norm()  = " << v2_new->norm() << endl;
    cout  << " v2_tens_keep->norm() = " << v2_tens_keep->norm() << endl << endl;
    
    v2_tens = v2_new;
#endif

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
