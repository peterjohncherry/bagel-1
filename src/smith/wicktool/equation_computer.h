#ifndef __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H
#define __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic.h>
#include <src/smith/wicktool/states_info.h>


namespace bagel {
namespace SMITH { 

namespace Equation_Computer { 

class Equation_Computer { 

    public: 
    Equation_Computer( std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<double>> eqn_info_in,
                       std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> Data_map_in );
    ~Equation_Computer(){};
  
    int nelea_ ;
    int neleb_ ;
    int ncore_ ;
    int norb_  ;
    int nstate_;
    int maxtile;
    int cimaxblock;
    std::shared_ptr<IndexRange> virt_  ;
    std::shared_ptr<IndexRange> active_  ;
    std::shared_ptr<IndexRange> closed_  ;
    std::shared_ptr<std::map< int , std::shared_ptr<IndexRange>>> ci_idxrng_map;

    std::shared_ptr<const Dvec> cc_; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> GammaMap;


    std::shared_ptr<std::map< std::string, std::shared_ptr<Dvec>>> dvec_sigma_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Determinants>>> det_old_map  ;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Civec>>> cvec_old_map; 

    std::shared_ptr<Equation<double>> eqn_info;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;

    std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> Tensor_Calc;

    /////////// Tensor contraction routines /////////////////////////

    //only for two contracted indexes
    std::shared_ptr<Tensor_<double>> contract_on_same_tensor( std::string Tname, std::string Tout_name, std::pair<int,int> ctr_todo) ;

    //for an arbitrary number of contracted indexes
    std::shared_ptr<Tensor_<double>> contract_on_same_tensor( std::string Tens_in, std::shared_ptr<std::vector<int>> contracted_index_positions ) ;

    std::shared_ptr<Tensor_<double>> contract_different_tensors( std::string T1name, std::string T2name, std::string Tout_name, std::pair<int,int> ctr_todo );

    std::shared_ptr<Tensor_<double>> reorder_block_Tensor(std::string Tname, std::shared_ptr<std::vector<int>> new_order);

    std::shared_ptr<Tensor_<double>> get_block_Tensor(std::string Tname);
  
    std::shared_ptr<Tensor_<double>> get_uniform_Tensor(std::shared_ptr<std::vector<std::string>> unc_ranges, double XX );

    ////////////Gamma routines (RDM class based) //////////////////

    std::shared_ptr<std::vector<std::shared_ptr<Tensor_<double>>>>
    get_gamma( std::string gamma_name );
    
    void get_wfn_data( std::shared_ptr<CIVecInfo<double>>  cvec_info );


    void compute_sigma2( std::string II_name, std::string JJ_name, std::string sigma_name );

    void compute_sigmaN( std::string predecessor_name, std::string gamma_name );
  
    void get_gamma2_from_sigma2_and_civec( std::string IBra_name,  std::string JKet_name );

    void sigma_2a1( double* cvec_ptr, double* sigma_ptr, std::shared_ptr<Determinants> dets );

    void sigma_2a2( double* cvec_ptr, double* sigma_ptr, std::shared_ptr<Determinants> dets ) ;

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
