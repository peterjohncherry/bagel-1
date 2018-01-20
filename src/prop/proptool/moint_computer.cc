#include <bagel_config.h>
#include <src/prop/proptool/moint.h>
#include <src/smith/wicktool/tensor_arithmetic.h>
#include <src/prop/proptool/moint_computer.h>

using namespace std;
using namespace bagel; 
using namespace bagel::SMITH::Tensor_Arithmetic; 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2(const  vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_v2" << endl;
  //flipping of indexes due to conflicting order definitions with current moint routine
  vector<SMITH::IndexRange>  alt_ordered_blocks = { blocks[3], blocks[1], blocks[2], blocks[0] };

  K2ext<DataType> v2 = K2ext<DataType>( info_, coeffs_, alt_ordered_blocks );

  // again for flipping indexes
  shared_ptr<vector<int>> alt_to_norm_order =  make_shared<vector<int>>( vector<int>  { 3, 1, 2, 0 }) ;
  auto Tensor_Calc = make_shared<Tensor_Arithmetic<DataType>>();
  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic<DataType>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);

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
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] =  *(range_conversion_map_->at(alt_ordered_blocks[ii])); 
  K2ext<DataType> v2 =  K2ext<DataType>( info_, coeffs_, blocks );

  // again for flipping indexes
  shared_ptr<vector<int>> alt_to_norm_order =  make_shared<vector<int>>( vector<int>  { 3, 1, 2, 0 } );
//  shared_ptr<SMITH::Tensor_<DataType>> v2_tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( v2.tensor(), alt_to_norm_order);
  return v2.tensor();// v2_tens; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//is the core fock minus diagonal component from above
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_h1( const vector<SMITH::IndexRange>& blocks, bool set_coeffs ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_h1" << endl;
//
  MOFock<DataType> h1( info_, blocks );

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
cout << "MOInt_Computer<DataType>::get_h1 string ver" << endl;

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] = *(range_conversion_map_->at(blocks_str[ii])); 
     
  MOFock<DataType> h1( info_, blocks );
  if ( set_coeffs ) 
    coeffs_ = h1.coeff();  

  return h1.tensor();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class MOInt_Computer<double>;
template class MOInt_Computer<std::complex<double>>;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
