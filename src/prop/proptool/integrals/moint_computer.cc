#include <bagel_config.h>
#include <src/prop/proptool/integrals/moint.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/integrals/moint_computer.h>

using namespace std;
using namespace bagel; 
using namespace bagel::Tensor_Arithmetic; 
#define __DEBUG_PROPTOOL_MOINT_COMPUTER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const  vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_v2" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //flipping of indexes due to conflicting order definitions with current moint routine
//  vector<SMITH::IndexRange>  alt_ordered_blocks = { blocks[3], blocks[1], blocks[2], blocks[0] };
  vector<SMITH::IndexRange>  alt_ordered_blocks = blocks;

  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, alt_ordered_blocks );

  // again for flipping indexes
  vector<int> alt_to_norm_order = { 3, 1, 2, 0 };
  auto Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);

  //TODO Why is this like this!?! Shouldn't you return v2_tens? Or just have the function in the return statement?
  return v2.tensor();
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
  //vector<string>  alt_ordered_blocks = {  blocks_str[3], blocks_str[1], blocks_str[2] , blocks_str[0] } ;
  vector<string>  alt_ordered_blocks = {  blocks_str[0], blocks_str[1], blocks_str[2] , blocks_str[3] } ;

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
    blocks[ii] =  *(range_conversion_map_->at(alt_ordered_blocks[ii])); 


  
  cout << endl;
  MOInt::K2ext_new<DataType> v2 =  MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );

  { //TEST  
    vector<SMITH::IndexRange> act4_ranges(blocks_str.size());
    for ( int ii = 0 ; ii != blocks_str.size(); ii++ )
      act4_ranges[ii] =  *(range_conversion_map_->at("a")); 
   
    cout << "v2_->norm() = " << v2.tensor()->norm(); cout.flush();  cout << " v2->size() = " << v2.tensor()->size_alloc() << endl;
    auto v2_act = Tensor_Arithmetic_Utils::get_sub_tensor( v2.tensor(), act4_ranges ); 
    cout << "v2_act->norm() = " << v2_act->norm(); cout.flush();  cout << " v2_act->size() = " << v2_act->size_alloc() << endl;
  
    Tensor_Arithmetic_Utils::Print_Tensor( v2_act, "v2_act" );  
  } //END TEST

  // again for flipping indexes
  vector<int> alt_to_norm_order =  { 3, 1, 2, 0 };
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);

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
  MOInt::MOFock_new<DataType> h1( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = h1.coeff();  

  return h1.tensor();
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
     
  MOInt::MOFock_new<DataType> h1( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = h1.coeff();  

  return h1.tensor();
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
  
  MOInt::MOFock_new<DataType> fock( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = fock.coeff();  

  return fock.tensor();
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
  MOInt::MOFock_new<DataType> fock( info_, blocks );

  if ( set_coeffs ) 
    coeffs_ = fock.coeff();  

  return fock.tensor();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_test_tensor( const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_MOINT_COMPUTER
cout << "MOInt_Computer<DataType>::get_test_tensor string ver" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<vector<SMITH::IndexRange>> blocks = make_shared<vector<SMITH::IndexRange>>(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks->at(ii) = *(range_conversion_map_->at(blocks_str[ii])); 

  return  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_Tensor_row_major( blocks );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class MOInt_Computer<double>;
template class MOInt_Computer<std::complex<double>>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
