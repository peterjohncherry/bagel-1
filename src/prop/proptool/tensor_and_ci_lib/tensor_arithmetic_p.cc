#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/debugging_utils.h>

//#define __DEBUG_PROPTOOL_TENSOR_ARITHMETIC_VERBOSE
#define __DEBUG_PROPTOOL_TENSOR_ARITHMETIC
#ifdef  __DEBUG_PROPTOOL_TENSOR_ARITHMETIC
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_debug.h>
#endif

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::Tensor_Sorter;
using namespace bagel::Tensor_Arithmetic_P_Utils; 
using namespace WickUtils;
using namespace Debugging_Utils;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Example Input : ( A_ijkl , B_wxyz, [ [ 0, 1, 2, 3] , [ 3, 2, 1, 0] ] , [ 1.0 , 2.0 ]   ) 
//Then will do A_ijkl += (1.0 * B_wxyz) 
//             A_ijkl += (2.0 * B_zyxw)
//The target data is changed, the summand is unchanged
//The inputs A and B must have equivalent ranges up to an reordering (e.g., A_[r1,r2] , B_r2,r1 is OK, A_[r1,r2,] 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void Tensor_Arithmetic_P::Tensor_Arithmetic_P<DataType>::add_list_of_reordered_tensors( shared_ptr<Tensor_<DataType>>& target,
                                                                                    shared_ptr<Tensor_<DataType>>& summand,
                                                                                    vector<vector<int>>& summand_reorderings,
                                                                                    vector<DataType>& summand_factors                      ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_TENSOR_ARITHMETIC 
cout << endl <<  "Tensor_Arithmetic_P::add_list_of_reordered_tensors" <<endl;  
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> target_ranges  = target->indexrange();
  vector<IndexRange> summand_ranges = summand->indexrange();

  cout.precision(13);

  if ( target->size_alloc() != 1 ) {
    vector<int> summand_maxs = get_num_index_blocks_vec( summand_ranges );
    vector<int> summand_mins(summand_maxs.size(), 0 );
    vector<int> summand_block_pos(summand_maxs.size(), 0 );
   
    do {

      vector<Index> summand_block_ranges = get_rng_blocks( summand_block_pos, summand_ranges );
      vector<Index> target_block_ranges = reorder_vector( summand_reorderings.front(), summand_block_ranges ) ;
       
      if ( target->is_local( target_block_ranges ) { 

        unique_ptr<DataType[]> target_block_data  = target->get_block(target_block_ranges);
        size_t summand_block_size  = summand->get_size( summand_block_ranges );    

        unique_ptr<DataType[]> summand_block_data  = summand->get_block(summand_block_ranges);
        DataType* summand_block_ptr = summand_block_data.get();
        vector<vector<int>>::iterator sr_it = summand_reorderings.begin();
        for ( typename vector<DataType>::iterator sf_it = summand_factors.begin(); sf_it !=  summand_factors.end() ; sr_it++, sf_it++ ) 
          unique_ptr<DataType[]> summand_block_data_reordered = reorder_tensor_data( summand_block_ptr, *sr_it, summand_block_ranges );

        target->put_block( target_block_data, target_block_ranges );
      }


    } while( fvec_cycle_skipper_f2b(summand_block_pos, summand_maxs, summand_mins) );

  } else { 
     for ( typename vector<DataType>::iterator sf_it = summand_factors.begin(); sf_it !=  summand_factors.end() ; ++sf_it ) 
       target->ax_plus_y( *sf_it , summand ); 
  } 

  return;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Tensor_Arithmetic_P::Tensor_Arithmetic_P<double>;
template class Tensor_Arithmetic_P::Tensor_Arithmetic_P<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
