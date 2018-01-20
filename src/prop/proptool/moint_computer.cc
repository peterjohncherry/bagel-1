#include <bagel_config.h>
#include <src/prop/proptool/moint.h>
#include <src/prop/proptool/moint_computer.h>

using namespace std;
using namespace bagel; 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//note, this does not have the diagonal component
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2( const vector<SMITH::IndexRange>& blocks ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_v2" << endl;
  K2ext<DataType> v2_ =  K2ext<DataType>( info_, coeffs_, blocks );
  return v2_.tensor(); ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Dupe routine with string input
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<SMITH::Tensor_<DataType>> MOInt_Computer<DataType>::get_v2( const vector<string>& blocks_str ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MOInt_Computer<DataType>::get_v2 string ver" << endl;

  vector<SMITH::IndexRange> blocks(blocks_str.size());
  for ( int ii = 0 ; ii != blocks_str.size(); ii++ ) 
    blocks[ii] =  *(range_conversion_map_->at(blocks_str[ii])); 
  K2ext<DataType> v2_ =  K2ext<DataType>( info_, coeffs_, blocks );
  return  v2_.tensor(); 
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
