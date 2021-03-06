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
  
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = calculate_v2_smith() ;
  v2_from_smith_ = v2_tens->copy();
  vector<int> alt_to_norm_order = { 3, 1, 2, 0};
  Debugging_Utils::print_sizes( v2_tens->indexrange(), "v2_range sizes" ); cout << endl;
  v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2_tens, alt_to_norm_order );

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

  MOInt::MOFock_new<DataType> one_el_ints( moint_info_, blocks );
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
cout << "MOInt_Computer<DataType>::calculate_fock : set_coeffs = " << set_coeffs << endl ;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  MOInt::MOFock_new<DataType> one_el_ints( moint_info_, blocks );
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

  bool  use_smith = false;
  if ( !use_smith ) { 
 
    shared_ptr<map<string,shared_ptr<SMITH::IndexRange>>>& rcm = range_conversion_map_; 
    auto rng = [&rcm] ( string range_name ) { return *(rcm->at(range_name)); };

    
    if ( !got_fock_coeffs_ )
      calculate_fock( {rng("free"), rng("free")}, true, true);
   
    //TODO excessive, for testing 
    vector<SMITH::IndexRange>  blocks = { rng("free"), rng("free"), rng("free"), rng("free") };
    MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( moint_info_, coeffs_, blocks );
    auto v2_tens = v2.tensor();
    v2_tens->allocate();
    v2_tens->zero();
    
    {
      vector<SMITH::IndexRange> keep_ranges_norder = { rng("free"), rng("free"), rng("free"), rng("free") };
      vector<SMITH::IndexRange> keep_ranges_sorder = { keep_ranges_norder[3], keep_ranges_norder[1], keep_ranges_norder[2], keep_ranges_norder[0] };
      shared_ptr<SMITH::Tensor_<DataType>> v2_keep =  v2.get_v2_part( keep_ranges_sorder ); 
      Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( v2_keep, v2_tens ); 
    }
    
    //TEST 
    v2_smith_ = v2_tens->copy(); 
    //END TEST
    return v2_tens;

  // TODO check this reordering!!
  } else {
    vector<int> r1032 = { 1, 0, 3, 2 };
    vector<int> r3120 = { 3, 1, 2, 0 };
    vector<int> reordering = r3120;
    auto v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2_from_smith_, reordering );
    Debugging_Utils::print_sizes( v2_from_smith_->indexrange(), "v2_from_smith rngs "); cout << endl;
    Debugging_Utils::print_sizes( v2_->indexrange(), "v2_ rngs "); cout << endl;
    return v2_; 
  }
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
cout << endl << "MOInt_Computer<DataType>::build_s_test_tensor "; cout.flush(); 
WickUtils::print_vector( ordering, "s_test_tensor_ordering"); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  bool use_smith = true;
  if ( use_smith ) { 
    Debugging_Utils::print_sizes( t_from_smith_->indexrange(), " t_from_smith" ); cout << endl;
    vector<int> r3120 = { 3, 1, 2, 0 };
    return Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( t_from_smith_, r3120 );

  } else { 
  
  shared_ptr<map<string,shared_ptr<SMITH::IndexRange>>>& rcm = range_conversion_map_; 
  auto rng = [&rcm] ( string range_name ) { return *(rcm->at(range_name)); };

  // Smith order : {3, 1, 2, 0}  
  vector<SMITH::IndexRange> full_block = { rng("notcore"), rng("notcore"), rng("notvirt"), rng("notvirt") };
  shared_ptr<SMITH::Tensor_<DataType>> s_test_tensor =  make_shared<SMITH::Tensor_<DataType>>( full_block ); 
  s_test_tensor->allocate();
  s_test_tensor->zero();

  return s_test_tensor; 
  }
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
