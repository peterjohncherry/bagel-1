#ifndef __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H
#define __SRC_SMITH_WICKTOOL_EQN_COMPUTER_H

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/util/kramers.h>
#include <src/smith/storage.h>
#include <src/smith/wicktool/tensor_sorter.h>

namespace bagel {
namespace SMITH { 

namespace Equation_Computer { 

class Equation_Computer { 

    public: 
    Equation_Computer(std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<double>> eqn_info,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>> >> CTP_data_map,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map);
    ~Equation_Computer(){};
  
    int nelea_ ;
    int neleb_ ;
    int ncore_ ;
    int norb_  ;
    int nstate_;
    std::shared_ptr<IndexRange> virt_  ;
    std::shared_ptr<IndexRange> active_  ;
    std::shared_ptr<IndexRange> closed_  ;
    std::shared_ptr<const Dvec> cc_; 
    std::shared_ptr<const Determinants> det_ ; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map;
    std::shared_ptr<std::unordered_map< std::string, std::shared_ptr<GammaInfo>>> GammaMap;

    std::shared_ptr<Equation<double>> eqn_info;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;

    ////////////////////////////////////
    std::shared_ptr<Tensor_<double>>
    get_block_Tensor(std::string Tname);

    std::shared_ptr<Tensor_<double>>
    contract_on_same_tensor( std::pair<int,int> ctr_todo, std::string Tname, std::string Tout_name) ;

    std::shared_ptr<Tensor_<double>>
    contract_different_tensors(std::pair<int,int> ctr_todo, std::string T1name, std::string T2name, std::string Tout_name);

    template<class DataType, class DType>
    std::shared_ptr<DType> contract_different_tensors( std::string T1name, std::string T2name,  std::pair<int,int> ctr_todo,
                                                       std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap ) ;
    //////////////////////////////
    std::shared_ptr<std::vector<int>> get_Tens_strides ( std::vector<int>& range_sizes) ;
  
    std::shared_ptr<std::vector<int>> get_CTens_strides( std::shared_ptr<std::vector<int>> range_sizes, int ctr1 , int ctr2 ) ;

    std::shared_ptr<std::vector<int>> get_CTens_strides( std::vector<int>& range_sizes, int ctr1 , int ctr2 ) ;

    std::shared_ptr<Tensor_<double>> get_uniform_Tensor(std::shared_ptr<std::vector<std::string>> unc_ranges, double XX );

    void Print_Tensor(std::shared_ptr<Tensor_<double>> Tens) ;

    std::vector<int> get_sizes(std::vector<Index>& Idvec) ;

    void Calculate_CTP(std::string A_contrib_name );

    template<class vtype>
    std::shared_ptr<std::vector<vtype>>
    inverse_reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) ;

    template<class vtype>
    std::shared_ptr<std::vector<vtype>>
    reorder_vector(std::shared_ptr<std::vector<int>> neworder , std::shared_ptr<std::vector<vtype>> origvec ) ;

    std::shared_ptr<std::vector<IndexRange>>
    Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    
    std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>>
    Get_Bagel_const_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str, std::shared_ptr<std::vector<int>> unc_pos);

    std::tuple< size_t, size_t >
    get_block_info(std::shared_ptr<std::vector<IndexRange>> id_ranges, std::shared_ptr<std::vector<int>> block_pos) ;

    void
    build_index_conversion_map(std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<IndexRange>>>> range_conversion_pairs );

    std::shared_ptr<std::vector<Index>>
    get_rng_blocks(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> old_ids) ;
   
    std::shared_ptr<std::vector<Index>>
    get_rng_blocks(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr< IndexRange>>> old_ids) ;

    std::vector<Index>
    get_rng_blocks(std::shared_ptr<std::vector<int>> forvec, std::vector<IndexRange>& id_ranges) ;
  
    std::vector<Index>
    get_rng_blocks_raw(std::shared_ptr<std::vector<int>> forvec, std::shared_ptr<std::vector<std::shared_ptr< IndexRange>>> old_ids) ;

    std::shared_ptr<std::vector<size_t>>
    get_sizes(std::shared_ptr<std::vector<Index>> Idvec);

    std::shared_ptr<std::vector<int>>
    get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec) ;

    std::shared_ptr<std::vector<size_t>>
    get_sizes(std::shared_ptr<std::vector<Index>> Idvec, int skip_id);

    std::shared_ptr<std::vector<int>>
    get_sizes(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec, int skip_id) ;

    size_t get_block_size(std::vector<Index>::iterator beginpos, std::vector<Index>::iterator endpos  ); 

    size_t get_unc_block_size( std::vector<Index>& idvec, std::pair<int,int> ctr ) ;

    size_t get_block_size(std::shared_ptr<std::vector<Index>> Idvec, int startpos, int endpos) ;

    std::shared_ptr<std::vector<int>>
    put_ctr_at_front(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);
    
    std::shared_ptr<std::vector<int>>
    put_ctr_at_back(std::shared_ptr<std::vector<int>> orig_pos , int ctr_pos);
    
    std::shared_ptr<Tensor_<double>>
    find_or_get_CTP_data(std::string CTP_name);

    std::shared_ptr<std::vector<int>>
    get_num_index_blocks_vec(std::shared_ptr<std::vector<std::shared_ptr<const IndexRange>>> rngvec) ;
 
    std::vector<int>
    get_num_index_blocks_vec(std::vector<IndexRange>& rngvec);
   
    std::pair<int,int>
    relativize_ctr_positions(std::pair <int,int> ctr_todo, std::shared_ptr<CtrTensorPart<double>>  CTP1,
                                                           std::shared_ptr<CtrTensorPart<double>>  CTP2);
    template<class DataType>
    std::unique_ptr<DataType[]>
    reorder_tensor_data( const DataType* orig_data, std::shared_ptr<std::vector<int>>  new_order_vec,
                         std::shared_ptr<std::vector<Index>> orig_index_blocks ) ;
    
    template<int N> 
    std::pair<std::vector<int>, std::pair<double,bool>> find_permutation(const KTag<N>& tag) const ;

    std::unique_ptr<double[]> get_block_of_data( double* data_ptr, std::shared_ptr<std::vector<IndexRange>> id_ranges, 
                                                              std::shared_ptr<std::vector<int>> block_pos) ;

    std::shared_ptr<std::vector<std::shared_ptr<Tensor_<double>>>>
    get_gammas(int MM , int NN, std::string gamma_name);
    
    std::shared_ptr<std::vector<std::shared_ptr<VectorB>>> compute_gammas(const int MM, const int NN ) ;

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
