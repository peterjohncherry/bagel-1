#ifndef __SRC_PROP_PROPTOOL_TENSOPCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TENSOPCOMPUTER_H

#include <tuple>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/task_translator/tensor_algop_info.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>
#include <src/prop/proptool/integrals/moint_computer.h>

namespace bagel {

namespace TensOp_Computer { 

template<class DataType>
class TensOp_Computer { 

    public: 
    TensOp_Computer( std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >> ACompute_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart_Base>>> CTP_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map,
                     std::shared_ptr<MOInt_Computer<DataType>> moint_computer ):
                     ACompute_map_(ACompute_map), CTP_map_(CTP_map), tensop_data_map_(tensop_data_map),
                     range_conversion_map_(range_conversion_map), moint_computer_(moint_computer), 
                     Tensor_Calc_(std::make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>()){}

    ~TensOp_Computer(){};
  
    std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >> ACompute_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart_Base>>> CTP_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> tensop_data_map_;

    std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<DataType>> Tensor_Calc_;
    std::shared_ptr<MOInt_Computer<DataType>> moint_computer_;

    /////////// Tensor contraction routines /////////////////////////

    void add_tensors( std::string target_name, std::string summand_name, DataType factor );
 
    void build_tensor( std::string new_data_name, std::vector<std::string> id_ranges, DataType init_val = (DataType)(0.0) ); 

    std::vector<SMITH::IndexRange> Get_Bagel_IndexRanges( const std::vector<std::string>& ranges_str ); 

    void sum_different_orderings( std::string target_name, std::string summand_name, std::vector<std::pair<double,double>> factor_list, std::vector<std::vector<int>> id_orders );

    //only for two contracted indexes
    std::shared_ptr<SMITH::Tensor_<DataType>> contract_on_same_tensor( std::string Tname, std::string Tout_name, std::pair<int,int> ctr_todo ) ;

    //for an arbitrary number of contracted indexes
    std::shared_ptr<SMITH::Tensor_<DataType>> contract_on_same_tensor( std::string Tens_in, std::shared_ptr<std::vector<int>> contracted_index_positions ) ;

    std::shared_ptr<SMITH::Tensor_<DataType>> contract_different_tensors( std::string T1name, std::string T2name, std::string Tout_name, std::pair<int,int> ctr_todo );

    std::shared_ptr<SMITH::Tensor_<DataType>> contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                                 std::pair<std::vector<int>,std::vector<int>> ctrs_todo                  );

    std::shared_ptr<SMITH::Tensor_<DataType>> direct_product_tensors( std::vector<std::string>& Tensor_names );

    std::shared_ptr<SMITH::Tensor_<DataType>> reorder_block_Tensor(std::string Tname, std::vector<int>& new_order);

    std::shared_ptr<SMITH::Tensor_<DataType>> get_block_Tensor(std::string Tname);

    void get_tensor_data_blocks(std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks );
  
    std::shared_ptr<SMITH::Tensor_<DataType>> get_uniform_Tensor(std::shared_ptr<std::vector<std::string>> unc_ranges, DataType XX );

    std::shared_ptr<SMITH::Tensor_<DataType>> divide_tensors( std::string T1_name, std::string T2_name );

    void divide_tensors_in_place( std::string T1_name, std::string T2_name );


    /////////// Utility routines /////////////////////////

    void Calculate_CTP( AContribInfo_Base& A_contrib_name );

    std::shared_ptr<std::vector<SMITH::IndexRange>>
    Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const SMITH::IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str, std::shared_ptr<std::vector<int>> unc_pos);

    void
    build_index_conversion_map(std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<SMITH::IndexRange>>>> range_conversion_pairs );

    std::shared_ptr<SMITH::Tensor_<DataType>> 
    find_or_get_CTP_data(std::string CTP_name);

    std::pair<int,int>
    relativize_ctr_positions(std::pair <int,int> ctr_todo, std::shared_ptr<CtrTensorPart_Base>  CTP1,
                                                           std::shared_ptr<CtrTensorPart_Base>  CTP2);

    void build_mo_tensor( std::string mo_tensor_name );

    void get_sub_tensor( std::string full_tens_name, std::string block_name,  const std::vector<std::string>& range_names );

    void build_test_tensor( std::string test_tensor_name, std::vector<size_t> dimensions, std::vector<size_t> max_blocks_sizes );

};

}
}
#endif
