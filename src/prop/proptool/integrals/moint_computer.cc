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
  v2_ = v2_tens;
//  vector<int> alt_to_norm_order = { 3, 1, 2, 0 };
//  v2_ = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2_tens, alt_to_norm_order );

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
 
  SMITH::IndexRange act = *(range_conversion_map_->at("a")); 
  SMITH::IndexRange core = *(range_conversion_map_->at("c")); 
  SMITH::IndexRange virt = *(range_conversion_map_->at("v")); 
  SMITH::IndexRange not_core = *(range_conversion_map_->at("notcor")); 
  SMITH::IndexRange not_virt = *(range_conversion_map_->at("notvir")); 
  SMITH::IndexRange free = core; free.merge(act); free.merge(virt);
   
  cout << "v2 ranges : core.size() = " << core.size() ; cout << " act.size() = " << act.size(); cout << " virt.size() = " << virt.size() << endl;

  if ( !got_fock_coeffs_ )
    calculate_fock( {free, free}, true, true);

  vector<SMITH::IndexRange>  blocks = { free, free, free, free };
  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );
  auto v2_tens = v2.tensor();
  v2_tens->allocate();
  v2_tens->zero();

  {
    vector<SMITH::IndexRange> keep_ranges_norder = { virt, virt, core, act };
    vector<SMITH::IndexRange> keep_ranges_sorder = { keep_ranges_norder[3], keep_ranges_norder[1], keep_ranges_norder[2], keep_ranges_norder[0] };
    shared_ptr<SMITH::Tensor_<DataType>> v2_keep =  v2.get_v2_part( keep_ranges_sorder ); 
    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( v2_keep, v2_tens ); 
  }
  {
    vector<SMITH::IndexRange> keep_ranges_norder = { virt, virt, act, core };
    vector<SMITH::IndexRange> keep_ranges_sorder = { keep_ranges_norder[3], keep_ranges_norder[1], keep_ranges_norder[2], keep_ranges_norder[0] };
    shared_ptr<SMITH::Tensor_<DataType>> v2_keep =  v2.get_v2_part( keep_ranges_sorder ); 
    Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( v2_keep, v2_tens ); 
  }

  cout << setprecision(13) <<  "v2_tens->norm() = "; cout.flush(); cout << v2_tens->norm() << endl;
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
cout << endl << "MOInt_Computer<DataType>::build_s_test_tensor "; cout.flush(); 
WickUtils::print_vector( ordering, "s_test_tensor_ordering"); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  SMITH::IndexRange act = *(range_conversion_map_->at("a")); 
  SMITH::IndexRange core = *(range_conversion_map_->at("c")); 
  SMITH::IndexRange virt = *(range_conversion_map_->at("v")); 

  // vector<int> alt_to_norm_order = { 3, 1, 2, 0 };
  SMITH::IndexRange not_active = *(range_conversion_map_->at("notact"));
  SMITH::IndexRange not_core =  *(range_conversion_map_->at("notcor"));
  SMITH::IndexRange not_virt =  *(range_conversion_map_->at("notvir"));
  SMITH::IndexRange free =  *(range_conversion_map_->at("free"));

//  Smith order : {3, 1, 2, 0}  
  vector<SMITH::IndexRange> full_block = { not_virt, not_virt, not_core, not_core };
  Debugging_Utils::print_sizes( WickUtils::reorder_vector( ordering, full_block ), "WickUtils::reorder_vector( ordering, full_block ) " ); cout << endl;
  shared_ptr<SMITH::Tensor_<DataType>> s_test_tensor =  make_shared<SMITH::Tensor_<DataType>>( WickUtils::reorder_vector( ordering, full_block ) ); 
  s_test_tensor->allocate();
  s_test_tensor->zero();

  Debugging_Utils::print_sizes( full_block,                  " full_block                "); cout << endl;
  Debugging_Utils::print_sizes( s_test_tensor->indexrange(), " s_test_tensor->id_ranges()"); cout << endl;

  {
  vector<SMITH::IndexRange> non_zero_block = { act, core, virt, virt }; 
  Debugging_Utils::print_sizes( WickUtils::reorder_vector( ordering, non_zero_block ), "WickUtils::reorder_vector( ordering, non_zero_block ) " ); cout << endl;
  shared_ptr<SMITH::Tensor_<DataType>> s_tmp =  make_shared<SMITH::Tensor_<DataType>>( WickUtils::reorder_vector( ordering, non_zero_block ) ); 
  s_tmp->allocate();
  s_tmp->zero();
  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_tmp, (DataType)(1.0) );
//  std::shared_ptr<SMITH::Tensor_<DataType>> s_tmp_as = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_tensor_antisymmetric( non_zero_block, (DataType)(1.0));
//  std::shared_ptr<SMITH::Tensor_<DataType>> s_tmp = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor(s_tmp_as, ordering );
  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( s_tmp, s_test_tensor ); 
  cout << " s_test_tensor->norm() * s_test_tensor->norm() = " <<  s_test_tensor->norm() * s_test_tensor->norm() <<endl;
  }
 {
  vector<SMITH::IndexRange> non_zero_block = { core, act, virt, virt }; 
  Debugging_Utils::print_sizes( WickUtils::reorder_vector( ordering, non_zero_block ), "WickUtils::reorder_vector( ordering, non_zero_block ) " ); cout << endl;
  shared_ptr<SMITH::Tensor_<DataType>> s_tmp =  make_shared<SMITH::Tensor_<DataType>>( WickUtils::reorder_vector( ordering, non_zero_block ) ); 
  s_tmp->allocate();
  s_tmp->zero();
  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_tmp, (DataType)(1.0) );
 // std::shared_ptr<SMITH::Tensor_<DataType>> s_tmp_as = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_uniform_tensor_antisymmetric( non_zero_block, (DataType)(-1.0));
 // std::shared_ptr<SMITH::Tensor_<DataType>> s_tmp = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor(s_tmp_as, ordering );
  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::put_sub_tensor( s_tmp, s_test_tensor ); 
  cout << " s_test_tensor->norm() * s_test_tensor->norm() = " <<  s_test_tensor->norm() * s_test_tensor->norm() <<endl;
  }

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
//  {
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, core, active }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(1.0) );
//  }{
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, active, core }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(-1.0) );
//  }{
//  vector<SMITH::IndexRange> non_zero_block = { virt, virt, active, active }; 
//  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( s_test_tensor, WickUtils::reorder_vector( ordering, non_zero_block ), (DataType)(0.0) );
//
// 
// 
// 
// 
// 
// 
// TESTING FOR BLOCK FETCHING
//
// vector<SMITH::IndexRange>  blocks = { free, free, free, free };
// MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );
// auto v2_tens = v2.tensor();
//
// vector<string> r1_names = { "c", "a", "notvir", "notcor", "v" };                                                                                              
// vector<string> r2_names = { "a", "v", "notcor", "notvir", "c" };                                                                                              
// vector<string> range_names(4);
// for ( auto r1_name : r1_names ) {                                                                                                                             
//   for ( auto r2_name : r2_names ){                                                                                                                            
//     range_names = { r1_name, r2_name, r1_name, r2_name } ;                                                                                                    
//     SMITH::IndexRange rng1 =  *( range_conversion_map_->at(r1_name) );                                                                                        
//     SMITH::IndexRange rng2 =  *( range_conversion_map_->at(r2_name) );                                                                                        
//     vector<SMITH::IndexRange> test_ranges = { rng1, rng2, rng1, rng2  };                                                                                      
//     WickUtils::print_vector( range_names, "test_range_names" ); Debugging_Utils::print_sizes( test_ranges, "test_range_sizes" ); cout << endl; 
//     auto v2_keep_orig = Tensor_Arithmetic_Utils::get_sub_tensor( v2_tens, test_ranges );                                                                      
//     shared_ptr<SMITH::Tensor_<DataType>> v2_keep =  v2.get_v2_part( test_ranges ); 
//     cout << " v2_keep_orig->norm() = "<<  v2_keep_orig->norm() ; cout.flush(); cout << " v2_keep->norm() = "<<  v2_keep->norm() <<endl;                       
//     v2_keep->ax_plus_y( -1.0, v2_keep_orig ); 
//     cout << " v2_keep_orig->norm() = "<<  v2_keep_orig->norm() ; cout.flush();  cout << " (v2_keep -  v2_keep_orig )->norm() = "<<  v2_keep->norm() << endl; 
//   }
// }
//
