#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_DEBUGGER_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_TENSOR_ARITH_DEBUGGER_H

#include <src/prop/proptool/proputils.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/debugging_utils.h>
#include <src/prop/proptool/proputils.h>
namespace bagel {

namespace Tensor_Arithmetic_Debugger { 

     template<typename DataType> 
     void add_tensor_along_trace_debug( std::shared_ptr<SMITH::Tensor_<DataType>> t_target, std::shared_ptr<SMITH::Tensor_<DataType>> t_summand,
                                        std::vector<int>& summand_pos, DataType factor ){ 
       std::cout << "Tensor_Arithmetic::add_tensor_along_trace_debug"; std::cout.flush();
       WickUtils::print_vector(summand_pos , "   summand_pos"); std::cout << std::endl;
       
       std::vector<SMITH::IndexRange> t_target_ranges_tmp  = t_target->indexrange();
       std::vector<SMITH::IndexRange> t_summand_ranges_tmp = t_summand->indexrange();
       {
       std::vector<int>::iterator sp_it = summand_pos.begin();
       for (std::vector<SMITH::IndexRange>::iterator tsr_it = t_summand_ranges_tmp.begin(); tsr_it != t_summand_ranges_tmp.end(); tsr_it++, sp_it++  )
         if  ( tsr_it->size() != t_target_ranges_tmp[*sp_it].size() ){
           std::cout << "t_summand range_block->size() = " << tsr_it->size() << " != " << t_target_ranges_tmp[*sp_it].size() << "  t_target_ranges[" << *sp_it <<"].size() "<< std::endl; 
           WickUtils::print_vector( summand_pos, "summand_pos") ; std::cout << std::endl;
           Debugging_Utils::print_sizes( t_target->indexrange(), "t_target_ranges_tmp" ); std::cout << std::endl;
           Debugging_Utils::print_sizes( t_summand->indexrange(), "t_summand_ranges_tmp" ); std::cout << std::endl;
           throw std::logic_error( "mismatched index range sizes in add_tensor_along_trace ! Aborting ! " ); 
         } 
       }
         
     };
  
     template<typename DataType> 
     void check_ranges( std::shared_ptr<SMITH::Tensor_<DataType>> t1, std::shared_ptr<SMITH::Tensor_<DataType>> t2 ) {
    
       std::vector<SMITH::IndexRange> t1_ranges = t1->indexrange();
       std::vector<SMITH::IndexRange> t2_ranges = t2->indexrange();
       
       for ( int ii = 0 ; ii != t2_ranges.size(); ii++ ) {  
         if ( t1_ranges.at(ii).size() != t2_ranges.at(ii).size() ) { 
           std::cout << "Ranges of t1 and t2 and not the same : " << std::endl; 
           Debugging_Utils::print_sizes( t1_ranges, "t1_range_sizes" ); std::cout << std::endl;
           Debugging_Utils::print_sizes( t2_ranges, "t2_range_sizes" ); std::cout << std::endl;
           throw std::logic_error( "Aborting!" ); 
         }
       }
      }


     template<class DataType>
     void check_all_same_ranges( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in ) { 
       if ( Tens_in->indexrange().size() == 0  )
         throw std::logic_error(" trying to take tensor of trace with no elements! Aborting!" );
 
       for ( int ii =1 ; ii != Tens_in->indexrange().size(); ii++ ) {
         if( Tens_in->indexrange()[ii].size() != Tens_in->indexrange()[ii-1].size() ){
           Debugging_Utils::print_sizes( Tens_in->indexrange(), "range_sizes" ); std::cout << std::endl;
           throw std::logic_error( " not all indexes of tensor have the same range, will not trace!! Aborting!" ); 
         }
       }
     }

     template<class DataType>
     void contract_on_same_tensor_debug( std::shared_ptr<SMITH::Tensor_<DataType>> Tens_in,  std::vector<int>& ctrs_pos) {
     std::cout << "Tensor_Arithmetic::contract_on_same_tensor" << std::endl;
     WickUtils::print_vector( ctrs_pos, "ctr_pos"); std::cout << std::endl; 
     Debugging_Utils::print_sizes( Tens_in->indexrange() , "Tens_in->indexrange()" ) ;std::cout << std::endl;
     { 
       std::vector<SMITH::IndexRange> id_ranges_in_tmp = Tens_in->indexrange();
       for ( int ii =1 ; ii != ctrs_pos.size(); ii++ ) {
         std::cout << id_ranges_in_tmp[ctrs_pos[ii]].size() << " " ; std::cout.flush();
         if( id_ranges_in_tmp[ctrs_pos[ii]].size() != id_ranges_in_tmp[ctrs_pos[ii-1]].size() ){
           std::cout << " trying to contract two index blocks of unequal lengths : " << std::endl;
           std::cout << " index at position  " <<  ctrs_pos[ii] << " has length " <<  id_ranges_in_tmp[ctrs_pos[ii]].size()  << std::endl; 
           std::cout << " index at position  " <<  ctrs_pos[ii-1] << " has length " <<  id_ranges_in_tmp[ctrs_pos[ii-1]].size() << std::endl;
           throw std::logic_error( "Aborting" ) ;
         } 
       } 
     }

    }
    template<typename DataType>
    void contract_different_tensors_debug( std::shared_ptr<SMITH::Tensor_<DataType>> Tens1_in, std::shared_ptr<SMITH::Tensor_<DataType>> Tens2_in,
                                           std::pair< std::vector<int>, std::vector<int> >& ctrs_todo ) {
       
      if (  ctrs_todo.first.size() !=  ctrs_todo.second.size() ) {
       std::cout << "  ctrs_todo.first.size() = " <<  ctrs_todo.first.size() << " != " <<  ctrs_todo.second.size() << " = ctrs_todo.second.size() " << std::endl;
       throw std::logic_error( "different number of contracted indexes on each tensor! Aborting! ");
      }
      std::vector<SMITH::IndexRange> T1_org_rngs_tmp = Tens1_in->indexrange();
      std::vector<SMITH::IndexRange> T2_org_rngs_tmp = Tens2_in->indexrange();
      for ( int ii = 0 ; ii != ctrs_todo.first.size(); ii++ ) {
        if ( T1_org_rngs_tmp[ctrs_todo.first[ii]].size() !=  T2_org_rngs_tmp[ctrs_todo.second[ii]].size() ){
          std::cout << " ctr1_size = " << T1_org_rngs_tmp[ctrs_todo.first[ii]].size();  std::cout.flush();
          std::cout << " ctr1_pos = " << ctrs_todo.first[ii]  << std::endl;
          std::cout << " ctr2_size = " << T2_org_rngs_tmp[ctrs_todo.second[ii]].size(); std::cout.flush();
          std::cout << " ctr2_pos = " << ctrs_todo.second[ii] << std::endl;
          throw std::logic_error("Extents of ranges to be contracted do not match!! Aborting");
        }
      }
    }
}
}
#endif
