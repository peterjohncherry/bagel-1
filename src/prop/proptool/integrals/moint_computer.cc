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
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_v2 IndexRange_ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //flipping of indexes due to conflicting order definitions with current moint routine
  vector<SMITH::IndexRange>  alt_ordered_blocks = blocks;

  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, alt_ordered_blocks );

  // again for flipping indexes
//  vector<int> alt_to_norm_order = { 3, 1, 2, 0 };
  vector<int> alt_to_norm_order = { 0, 1, 2, 3 };
  auto Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order );

  //TODO Why is this like this!?! Shouldn't you return v2_tens? Or just have the function in the return statement?
  cout << "MO_comp  v2_tens->norm() = " <<   v2_tens->norm() << endl;

  return v2_tens;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_v2 string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //flipping of indexes due to conflicting order definitions with current moint routine
  vector<string>  alt_ordered_blocks = {  blocks_str[0], blocks_str[1], blocks_str[2] , blocks_str[3] } ;

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
    blocks[ii] =  *(range_conversion_map_->at(alt_ordered_blocks[ii])); 
  
  MOInt::K2ext_new<DataType> v2 =  MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );

  // again for flipping indexes
//  vector<int> alt_to_norm_order =  { 3, 1, 2, 0 };
  vector<int> alt_to_norm_order =  { 0, 1, 2, 3 };
//  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = v2.tensor();
  {
    vector<SMITH::IndexRange> blocks_act(blocks_str.size());
    string act_name = "a";
    for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
      blocks_act[ii] =  *(range_conversion_map_->at(act_name)); 

    shared_ptr<SMITH::Tensor_<DataType>> v2_act = 
    Tensor_Arithmetic_Utils::get_sub_tensor( v2_tens, blocks_act );
    Tensor_Arithmetic_Utils::print_tensor_with_indexes( v2_act, " INT COMPUTER v2_act" );  cout << endl;
  }

  return  v2_tens; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is the core ham minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_h1( const vector<SMITH::IndexRange>& blocks, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_h1" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = one_el_ints.coeff();  

  return one_el_ints.h1();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_h1( const vector<string>& blocks_str, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_h1 string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 
     
  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = one_el_ints.coeff();  

  return one_el_ints.h1();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_fock( const vector<string>& blocks_str, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_fock string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 
  
  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = one_el_ints.coeff();  

  return one_el_ints.fock();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is core_fock minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_fock( const vector<SMITH::IndexRange>& blocks, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_fock" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  MOInt::MOFock_new<DataType> one_el_ints( info_, blocks );

  if ( set_coeffs ) 
    coeffs_ = one_el_ints.coeff();  

  return one_el_ints.fock();
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
