#include <bagel_config.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/tensor_and_ci_lib/vector_bundle.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;
//#define __DEBUG_PROPTOOL_VECTOR_BUNDLE

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
Vector_Bundle<DataType>::Vector_Bundle( vector<int> id_range_sizes, size_t vec_len, size_t vec_maxtile,
                                        bool all_sparse, bool alloc, bool zero ) {
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
Vector_Bundle<DataType>::set_vectors( DataType value ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::set_vectors" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   for ( auto elem : *vector_map_ )
     Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( elem.second , value);
   
   return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Vector_Bundle<DataType>::set_test_vectors() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::set_test_vectors" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   vector<int> powers_of_ten((vector_map_->begin()->first).size());
   powers_of_ten.back() =1;
   if ( powers_of_ten.size() >1 )
     for ( vector<int>::reverse_iterator pot_rit = powers_of_ten.rbegin()+1; pot_rit != powers_of_ten.rend(); pot_rit++ )
       *pot_rit = (*(pot_rit-1))*10;

   for ( auto& elem : *vector_map_ ){
     DataType value = (DataType)(inner_product( elem.first.begin(), elem.first.end(), powers_of_ten.begin(), 0));
     Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( elem.second , value);
   }  
   return;
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Vector_Bundle<DataType>::merge_fixed_ids( shared_ptr<Vector_Bundle<DataType>> bundle_to_merge, vector<int>& fixed_ids,
                                          vector<bool> ids_overwrite_pattern, char option, DataType factor ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//if fixed_ids = [1,2,3] , ids_overwrite_pattern = [ 0, 1, 1, 0, 0]
//then vector at [4, 9] in bundle_to_merge->vector_map_ is put at  [ 1, 4, 9, 2, 3] in current map 
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<DataType>::merge_fixed_ids" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( auto& new_bundle_elem : *(bundle_to_merge->vector_map()) ) {
    vector<int> new_ids(ids_overwrite_pattern.size());
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
    if ( new_ids_loc == vector_map_->end() ) {
      vector_map_->emplace(new_ids, new_bundle_elem.second );
    } else if (option == 'O') { 
      new_ids_loc->second = new_bundle_elem.second; 
    } else if (option == 'C') { 
      new_ids_loc->second->ax_plus_y( factor, new_bundle_elem.second ); 
    } else { 
      cout << "did not overwrite or combine; vector " ;cout.flush(); WickUtils::print_vector( new_ids ); cout << " is already in map" << endl; 
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void
Vector_Bundle<std::complex<double>>::print(string name, int line_max) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<std::complex<double>>::print" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   throw logic_error( " Cannot print Vector_Bundle<DataType> when DataType is complex<double>! Aborting; " );
   return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
void
Vector_Bundle<double>::print(string name ,  int line_max ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_VECTOR_BUNDLE
cout << "Vector_Bundle<double>::print" << endl;
//NOTE : This assumes all vectors have the same IndexRange, and will only print the first block.
//       Update at some point, but really should be fine for now
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (name != "" )
    cout << "---------------" << name <<" ------------------ "<< endl;

  vector<Index> block = { idx_range_.range(0) }; 
  auto block_size = block[0].size(); 
 
  for( auto vm_it = vector_map_->begin(); vm_it != vector_map_->end() ; vm_it++ ) {
    print_vector( vm_it->first );
    vector<SMITH::Index> range_blocks = vm_it->second->indexrange()[0].range();
    for ( const auto& block : range_blocks ){
      unique_ptr<double[]> data_block =  vm_it->second->get_block( vector<Index> { block } );
      double* end_ptr = data_block.get() + block.size();
      int ii = 0; 
      for ( double* data_ptr = data_block.get(); data_ptr != end_ptr; ++data_ptr )
        cout << *data_ptr << "  " ; 
      cout << endl;
      vm_it->second->put_block( data_block, vector<Index> {block} ); 
    }
  }

//TODO find out why the below code does not work; I suspect it is to do with fetching the pointers
//
//  auto vm_it = vector_map_->begin();
//  do {
//    vector<double*> ptr_vec(0);
//    int line_lim =0;
//    cout << endl << endl;
//    while ( vm_it != vector_map_->end() && line_lim != line_max ) {
    //  ptr_vec.push_back(vm_it->second->get_block( block ).get()); <---- THIS LINE IS PROBABLY THE PROBLEM
//      print_vector( vm_it->first );
//      ++line_lim;
//      ++vm_it;
//    }

//    cout << endl;
//    size_t elem_num = 0;
//    while ( elem_num != block_size ) {
//      for ( vector<double*>::iterator pv_it = ptr_vec.begin(); pv_it != ptr_vec.end(); ++pv_it ){ 
//        cout << *(*pv_it) << "   " ;
//        ++(*pv_it);
//      }
//      ++elem_num;
//      cout << endl;
//   }

//  } while ( vm_it != vector_map_->end() );

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Vector_Bundle<double>;
template class Vector_Bundle<std::complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
