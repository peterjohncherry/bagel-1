#ifndef __SRC_SMITH_WICKTOOL_INTERFACE_H
#define __SRC_SMITH_WICKTOOL_INTERFACE_H

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/smith.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>

namespace bagel {
namespace SMITH { 

namespace Equation_Computer { 

class Equation_Computer { 
    public: 
    Equation_Computer(std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<Tensor_<double>>> eqn_info,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>> >> CTP_data_map,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map);
    ~Equation_Computer(){};
  
    int  nelea_ ;
    int  neleb_ ;
    int  ncore_ ;
    int  norb_  ;
    int  nstate_;
    std::shared_ptr<IndexRange> virt_  ;
    std::shared_ptr<IndexRange> active_  ;
    std::shared_ptr<IndexRange> closed_  ;
    std::shared_ptr<const Dvec> cc_; 
    std::shared_ptr<const Determinants> det_ ; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<Tensor_<double>>>>> CTP_map;

    std::shared_ptr<Equation<Tensor_<double>>> eqn_info;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;

    template<class vtype>
    std::shared_ptr<std::vector<vtype>> inverse_reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) ;

    template<class vtype>
    std::shared_ptr<std::vector<vtype>> reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) ;


    std::shared_ptr<Tensor_<double>> contract_on_same_tensor( std::pair<int,int> ctr_todo, std::string Tname);

    std::shared_ptr<Tensor_<double>> contract_different_tensors(std::pair<int,int> ctr_todo, std::string T1name, std::string T2name );

    std::shared_ptr<std::vector<IndexRange>> Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);
  
    void build_index_conversion_map(std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<IndexRange>>>> range_conversion_pairs );

    template<class DataType, class DType>
    std::shared_ptr<DType> contract_different_tensors( std::string T1name, std::string T2name,  std::pair<int,int> ctr_todo,
                                                       std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap ) ;

    std::shared_ptr<std::vector<Index>> get_rng_blocks(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> old_ids) ;
    std::shared_ptr<std::vector<Index>> get_rng_blocks(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr< IndexRange>>> old_ids) ;
    std::vector<Index> get_rng_blocks_raw(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr< IndexRange>>> old_ids) ;

    std::shared_ptr<std::vector<size_t>> get_sizes(std::shared_ptr<std::vector<Index>> Idvec);

    std::shared_ptr<std::vector<int>> get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec) ;

    std::shared_ptr<Tensor_<double>> get_block_Tensor(std::string Tname);

    size_t get_block_size(std::shared_ptr<std::vector<Index>> Idvec, int startpos, int endpos) ;

    template<class DataType, class DType>
    std::unique_ptr<DataType[]>
    get_reordered_Tensor_data(std::shared_ptr<std::vector<int>> rng_block_pos, std::shared_ptr<std::vector<const IndexRange>> T_org_rng,
                              std::shared_ptr<std::vector<const IndexRange>> T_new_rng, std::shared_ptr<DType> Tens )  ;

    template<class DataType>
    std::unique_ptr<DataType[]>
    reorder_tensor_data(const DataType* orig_data,  size_t data_size, std::vector<int>  new_order_vec, std::vector<size_t> new_sizes_vec ) ;

    
    
    std::unique_ptr<double[]> get_block_of_data( double* data_ptr, std::shared_ptr<std::vector<IndexRange>> id_ranges, 
                                                              std::shared_ptr<std::vector<int>> block_pos) ;


    std::shared_ptr<std::vector<std::shared_ptr<Tensor_<double>>>>
    get_gammas(int MM , int NN, std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> gamma_ranges_REMOVE_AFTER_DEBUG);

    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>, std::shared_ptr<RDM<3>> >
    compute_gamma12(const int MM, const int NN ) ;
 
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>, std::shared_ptr<RDM<3>> >
    compute_gamma12_from_civec(std::shared_ptr<const Civec> cbra, std::shared_ptr<const Civec> cket) const ;
 
    std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>, std::shared_ptr<RDM<3>> >
    compute_gamma12_last_step(std::shared_ptr<const Dvec> dbra, std::shared_ptr<const Dvec> dket, std::shared_ptr<const Civec> cibra) const ;

    void sigma_2a1(std::shared_ptr<const Civec> cvec, std::shared_ptr<Dvec> sigma) const ;

    void sigma_2a2(std::shared_ptr<const Civec> cvec, std::shared_ptr<Dvec> sigma) const ;

};
}
}
}
////////////////////////////////////////////////////
//template class CtrTensorPart<std::vector<double>>;  
//template class CtrTensorPart<Tensor_<double>>;
///////////////////////////////////////////////////
#endif
