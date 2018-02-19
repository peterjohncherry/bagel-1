#include <bagel_config.h>
#include <src/prop/proptool/integrals/moint.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/integrals/moint_computer.h>

using namespace std;
using namespace bagel; 
using namespace bagel::Tensor_Arithmetic; 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const  vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_v2" << endl;
  //flipping of indexes due to conflicting order definitions with current moint routine
  vector<SMITH::IndexRange>  alt_ordered_blocks = { blocks[3], blocks[1], blocks[2], blocks[0] };

  MOInt::K2ext_new<DataType> v2 = MOInt::K2ext_new<DataType>( info_, coeffs_, alt_ordered_blocks );

  // again for flipping indexes
  shared_ptr<vector<int>> alt_to_norm_order =  make_shared<vector<int>>( vector<int>  { 3, 1, 2, 0 }) ;
  auto Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);

  return v2.tensor(); //  v2_tens;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_v2 string ver" << endl;

  cout << "blocks    = [ "; for ( string elem : blocks_str ) { cout << elem << " " ; } cout << "]" <<  endl;
  //flipping of indexes due to conflicting order definitions with current moint routine
  vector<string>  alt_ordered_blocks = {  blocks_str[3], blocks_str[1], blocks_str[2] , blocks_str[0] } ;
  cout << "blocks3120 = [ "; for ( string elem : alt_ordered_blocks ) { cout << elem << " " ; } cout << "]" <<  endl;

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ){ 
    blocks[ii] =  *(range_conversion_map_->at(alt_ordered_blocks[ii])); 
    cout << "blocks["<<ii<<"] = " <<  blocks[ii].size(); 
  }
  cout << endl;
  MOInt::K2ext_new<DataType> v2 =  MOInt::K2ext_new<DataType>( info_, coeffs_, blocks );

   cout << " v2.tensor->norm()->size_alloc() = " <<  v2.tensor()->size_alloc() << endl;
   cout << " v2.tensor()->norm() = " << v2.tensor()->norm() << endl; 
  // again for flipping indexes
  shared_ptr<vector<int>> alt_to_norm_order =  make_shared<vector<int>>( vector<int>  { 3, 1, 2, 0 } );
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);
  cout << "X" << endl;
  //return v2.tensor();// v2_tens; 
  return  v2_tens; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is the core ham minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_h1( const vector<SMITH::IndexRange>& blocks, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_h1" << endl;

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  cout << "pre fock coeffs_->norm() = " << coeffs_->norm() << endl;
  auto bob = make_shared<MatType> (*coeffs_) ;
  cout << "pre fock bob->norm() = " << bob->norm() << endl;
 
  MOInt::MOFock_new<DataType> h1( info_, blocks );
 
  bob->ax_plus_y( -1.0 , h1.coeff()); 

  cout << " coeffs_->ax_plus_y( -1.0 , h1.coeff())->norm() = " <<  bob->norm() << endl; 

  if ( set_coeffs ) 
    coeffs_ = h1.coeff();  

  cout << "new  coeffs_->norm() = " <<  coeffs_->norm() << endl; 
  return h1.tensor();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_h1( const vector<string>& blocks_str, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_h1 string ver" << endl;

  using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 
     
  cout << "pre fock coeffs_->norm() = " << coeffs_->norm() << endl;
  auto bob = make_shared<MatType> (*coeffs_) ;
  cout << "pre fock bob_->norm() = " << bob->norm() << endl;
  
  MOInt::MOFock_new<DataType> h1( info_, blocks );

  bob->ax_plus_y( -1.0 , h1.coeff()); 

  cout << " coeffs_->ax_plus_y( -1.0 , h1.coeff())->norm() = " <<  bob->norm() << endl; 

  if ( set_coeffs ) 
    coeffs_ = h1.coeff();  

  cout << "new  coeffs_->norm() = " <<  coeffs_->norm() << endl; 

  return h1.tensor();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_fock( const vector<string>& blocks_str, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_h1 string ver" << endl;

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
cout << "MOInt_Computer<DataType>::get_fock" << endl;

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
cout << "MOInt_Computer<DataType>::get_test_tensor string ver" << endl;

  shared_ptr<vector<SMITH::IndexRange>> blocks = make_shared<vector<SMITH::IndexRange>>(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks->at(ii) = *(range_conversion_map_->at(blocks_str[ii])); 

  return  Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_Tensor_row_major( blocks );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class MOInt_Computer<double>;
template class MOInt_Computer<std::complex<double>>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
