#ifndef __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H
#define __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H

#include <tuple>
#include <src/smith/wicktool/expression.h>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic.h>
#include <src/smith/wicktool/states_info.h>


namespace bagel {
namespace SMITH { 

namespace TensOp_Computer { 

class TensOp_Computer { 

    public: 
    TensOp_Computer( std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >> ACompute_map_in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Data_map_in );
    ~TensOp_Computer(){};
  
    std::shared_ptr<std::map< std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >> ACompute_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Data_map;



    std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> Tensor_Calc;

    /////////// Tensor contraction routines /////////////////////////

    //only for two contracted indexes
    std::shared_ptr<Tensor_<double>> contract_on_same_tensor( std::string Tname, std::string Tout_name, std::pair<int,int> ctr_todo ) ;

    //for an arbitrary number of contracted indexes
    std::shared_ptr<Tensor_<double>> contract_on_same_tensor( std::string Tens_in, std::shared_ptr<std::vector<int>> contracted_index_positions ) ;

    std::shared_ptr<Tensor_<double>> contract_different_tensors( std::string T1name, std::string T2name, std::string Tout_name, std::pair<int,int> ctr_todo );

    std::shared_ptr<Tensor_<double>> reorder_block_Tensor(std::string Tname, std::shared_ptr<std::vector<int>> new_order);

    std::shared_ptr<Tensor_<double>> get_block_Tensor(std::string Tname);
  
    std::shared_ptr<Tensor_<double>> get_uniform_Tensor(std::shared_ptr<std::vector<std::string>> unc_ranges, double XX );

    /////////// Utility routines /////////////////////////

    void Calculate_CTP(std::string A_contrib_name );

    std::shared_ptr<std::vector<IndexRange>>
    Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str, std::shared_ptr<std::vector<int>> unc_pos);

    void
    build_index_conversion_map(std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<IndexRange>>>> range_conversion_pairs );

    std::shared_ptr<Tensor_<double>>
    find_or_get_CTP_data(std::string CTP_name);

    std::pair<int,int>
    relativize_ctr_positions(std::pair <int,int> ctr_todo, std::shared_ptr<CtrTensorPart<double>>  CTP1,
                                                           std::shared_ptr<CtrTensorPart<double>>  CTP2);



};

}
}
}
#endif
