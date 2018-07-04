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
       cout << "all_sparse" << endl;
       ids = mins;
//     print_vector( maxs, "maxs"); print_vector( mins, "mins"); cout << endl;
       do { 
         sparsity_map_->emplace( ids, true );
//       print_vector(ids, "ids"); cout << endl;
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
   cout << "leaving init" << endl;
   return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Vector_Bundle<DataType>::merge_fixed_ids( shared_ptr<Vector_Bundle> bundle_to_merge, vector<int>& fixed_ids, vector<bool> ids_overwrite_pattern,
                                          bool overwrite ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//if fixed_ids = [1,2,3] , ids_overwrite_pattern = [ 0, 1, 1, 0, 0]
//then vector at [4, 9] in bundle_to_merge->vector_map_ is put at  [ 1, 4, 9, 2, 3] in current map 
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::merge_fixed_ids" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( auto& new_bundle_elem : *(bundle_to_merge->vector_map()) ) {
    vector<int> new_ids(new_bundle_elem.first.size());
    {
    vector<int>::iterator ni_it = new_ids.begin();
    vector<int>::iterator fi_it = fixed_ids.begin();
    vector<int>::const_iterator mi_it = new_bundle_elem.first.begin();
    for ( auto change_id : ids_overwrite_pattern ){
      if( !change_id ) { 
        *ni_it = *mi_it;
         ++mi_it;
      } else { 
        *ni_it = *fi_it;
        ++fi_it;
      }
      ++ni_it;
    }
    }
    auto new_ids_loc = vector_map_->find( new_ids );
    if ( new_ids_loc ==vector_map_->end() ) {
      vector_map_->emplace(new_ids, new_bundle_elem.second );
    } else if (overwrite) { 
      new_ids_loc->second = new_bundle_elem.second; 
    } else { 
      cout << "did not overwrite; vector " ;cout.flush(); WickUtils::print_vector( new_ids ); cout << " is already in map" << endl; 
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Vector_Bundle<double>;
template class Vector_Bundle<std::complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
