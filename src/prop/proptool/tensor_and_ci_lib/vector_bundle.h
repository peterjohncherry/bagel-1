#ifndef __SRC_PROP_PROPTOOL_VECTOR_BUNDLE
#define __SRC_PROP_PROPTOOL_VECTOR_BUNDLE

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <vector>
#include <map>
#include <memory>
#include <iostream>

namespace bagel {
// This is basically just like a stripped down Dvec, only each one of the civectors is in Tensor format.
// Also has sparsity info.
template<typename DataType>
class Vector_Bundle {  
 
    private: 
   
      std::shared_ptr< std::map< std::vector<int>, std::shared_ptr<SMITH::Tensor_<DataType>> >> vector_map_;
      std::shared_ptr< std::map< std::vector<int>, bool >> sparsity_map_;
      std::vector<int> id_range_sizes_;
      size_t vec_len_;
      size_t vec_maxtile_;
 
    public: 

      SMITH::IndexRange idx_range_;
      std::vector<SMITH::IndexRange> index_range_vec_;

      Vector_Bundle();
      Vector_Bundle(std::vector<int> id_range_sizes, size_t vec_len, size_t vec_maxtile,
                    bool all_sparse = true, bool alloc = false, bool zero = false ); 
      ~Vector_Bundle(){};

      void init( std::vector<int> id_range_sizes, size_t vec_len, size_t vec_maxtile, bool all_sparse, bool alloc, bool zero ); 

      void set_vector( std::vector<int>& ids , std::shared_ptr<SMITH::Tensor_<DataType>> vec, bool overwrite = true ) { 

        auto vmap_loc = vector_map_->find( ids ); 
              
        if ( vmap_loc != vector_map_->end() ) {
          if ( overwrite ) { 
            vmap_loc->second = vec;
            auto smap_loc = sparsity_map_->find( ids ); 
            if ( smap_loc != sparsity_map_->end() ) 
              smap_loc->second = false;

          } else { 
            std::cout << "vector at [ "; std::cout.flush();
            std::for_each( ids.begin(), ids.end(), [](int ii ) { std::cout << ii << " " ; } );
            std::cout.flush(); std::cout << " ] already in map" << std::endl;
          } 

        } else {
          vector_map_->emplace( ids, vec );
          auto smap_loc = sparsity_map_->find( ids ); 
          if ( smap_loc != sparsity_map_->end() ){ 
            smap_loc->second = false;
          } else {
            sparsity_map_->emplace(ids, false); 
          }
        }
      }

      void set_sparsity(std::vector<int>& ids, bool sparse = true ) {
        auto smap_loc = sparsity_map_->find( ids ); 
        if ( smap_loc != sparsity_map_->end() ){ 
          smap_loc->second = sparse;
        } else {
          sparsity_map_->emplace(ids, sparse); 
        }
        return;
      } 


      std::shared_ptr<SMITH::Tensor_<DataType>> vector_map(const std::vector<int>& ids ) const { return vector_map_->at(ids); } 

      bool sparsity_map(const std::vector<int>& ids ) const {
        auto smap_loc = sparsity_map_->find( ids ); 
        if ( smap_loc != sparsity_map_->end() ){ 
          return smap_loc->second;
        } else {
          return true;
        }
      }

      // TODO get rid of these
      std::shared_ptr< std::map< std::vector<int>, bool>> sparsity_map() { return sparsity_map_; }
      std::shared_ptr< std::map< std::vector<int>, std::shared_ptr<SMITH::Tensor_<DataType>> >> vector_map() { return vector_map_; }
};

}
#endif
