#include <bagel_config.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/tensor_and_ci_lib/vector_bundle.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;
#define __DEBUG_PROPTOOL_VECTOR_BUNDLE

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
Vector_Bundle<DataType>::Vector_Bundle() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::Vector_Bundle 1" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  sparsity_map_ = std::make_shared<std::map< std::vector<int>, bool >>();
  vector_map_ = std::make_shared<std::map< std::vector<int>, std::shared_ptr<SMITH::Tensor_<DataType>> >>();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
Vector_Bundle<DataType>::Vector_Bundle( vector<int> id_range_sizes, size_t vec_len, size_t vec_maxtile, bool all_sparse, bool alloc, bool zero ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::Vector_Bundle 2" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   sparsity_map_ = std::make_shared<std::map< std::vector<int>, bool >>();
   vector_map_ = std::make_shared<std::map< std::vector<int>, std::shared_ptr<SMITH::Tensor_<DataType>> >>();
   init( id_range_sizes, vec_len, vec_maxtile, all_sparse, alloc, zero );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Vector_Bundle<DataType>::init( vector<int> id_range_sizes, size_t vec_len, size_t vec_maxtile, bool all_sparse, bool alloc, bool zero ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::init" << endl;
if ( zero && !alloc ) 
  throw logic_error( "Cannot zero Vector_Bundle without allocating... Aborting! "); 
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   id_range_sizes_ = id_range_sizes; 
   vec_len_ = vec_len ;
   vec_maxtile_ = vec_maxtile;
   idx_range_ = IndexRange( vec_len_, vec_maxtile_ );
   index_range_vec_ = vector<IndexRange>( 1, idx_range_ );

   if ( all_sparse || alloc || zero ) {  

     vector<int> mins(id_range_sizes_.size(), 0);
     vector<int> maxs(id_range_sizes_.size());
     
     {
     vector<int>::iterator mx_it = maxs.begin();
     for ( vector<int>::iterator irs_it = id_range_sizes_.begin(); irs_it != id_range_sizes_.end(); ++irs_it, ++mx_it ) 
       *mx_it = (*irs_it)-1;
     
     }
     
     vector<int> ids;

     if ( all_sparse ) {
       ids = mins;
       do { 
         sparsity_map_->emplace( ids, true );
       } while ( fvec_cycle_skipper( ids, maxs, mins) );
     }
   
     if (alloc && !zero) {
       ids = mins;
       do { 
         shared_ptr<Tensor_<DataType>> new_vec = make_shared<Tensor_<DataType>>( index_range_vec_ ); 
         new_vec->allocate();
         vector_map_->emplace( ids, new_vec );
       } while ( fvec_cycle_skipper( ids, maxs, mins) );

     } else if ( alloc && zero ) { 

       ids = mins;
       do { 
         shared_ptr<Tensor_<DataType>> new_vec = make_shared<Tensor_<DataType>>( index_range_vec_ ); 
         new_vec->allocate();
         new_vec->zero();
         vector_map_->emplace( ids, new_vec );
       } while ( fvec_cycle_skipper( ids, maxs, mins) );
     }
   }
   return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Vector_Bundle<double>;
template class Vector_Bundle<std::complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
