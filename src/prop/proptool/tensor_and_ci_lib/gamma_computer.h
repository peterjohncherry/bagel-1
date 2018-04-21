#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_GAMMA_COMPUTER_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_GAMMA_COMPUTER_H

#include <tuple>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/ci/fci/civec.h>
#include <src/ci/fci/dvec.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
namespace bagel {

template< typename DataType> 
class Gamma_Computer { 

    public: 
      Gamma_Computer( std::shared_ptr< std::map< std::string, std::shared_ptr<GammaInfo<DataType>>>>          Gamma_info_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>    CIvec_data_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>    sigma_data_map_,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>    gamma_data_map_,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<const Determinants>>> Determinants_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<SMITH::IndexRange>>>         range_conversion_map ); 
      ~Gamma_Computer(){};
      
      int nelea_  ;
      int neleb_  ;
      int norb_   ;
      int nstate_ ;
      int maxtile ;
      int cimaxblock;
      std::shared_ptr<SMITH::IndexRange> virt_   ;
      std::shared_ptr<SMITH::IndexRange> active_ ;
      std::shared_ptr<SMITH::IndexRange> closed_ ;
      
      std::shared_ptr<const Dvec> cc_; 
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> CIvec_data_map;

      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo<DataType>>>> Gamma_info_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<const Determinants>>> Determinants_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<Determinants>>> Determinants_map_new;
      
      std::shared_ptr<Expression<DataType>> Eqn_info;
      
      std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<DataType>> Tensor_Calc;
      
      /////////// Gamma routines /////////////////////////

      void get_gamma_tensor( std::string gamma_name);
      
      void build_gamma2_tensor(std::shared_ptr<GammaInfo<DataType>> gamma2_info );

      void build_sigma2_tensor(std::shared_ptr<GammaInfo<DataType>> gamma2_info );

      void build_sigmaN_tensor(std::shared_ptr<GammaInfo<DataType>> gamma_info );

      void build_gammaN_tensor(std::shared_ptr<GammaInfo<DataType>> gamma_info );

      void build_sigmaN_block( std::string sigmaN,     std::vector<SMITH::Index>& id_blocks_ps, std::vector<int>& offsets_ps,
                               std::string prev_sigma, std::vector<SMITH::Index>& id_blocks_ij, std::vector<int>& offsets_ij  ) ;

      void sigma_2a1( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<const Determinants> dets, int cvec_offset, std::vector<int>& sigma_offsets, 
                      int cvec_block_size, std::vector<int>& sigma_block_sizes );

      void sigma_2a2( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<const Determinants> dets, int cvec_offset, std::vector<int>& sigma_offsets, 
                      int cvec_block_size, std::vector<int>& sigma_block_sizes );

      /////////// Utility routines /////////////////////////
      
      std::shared_ptr<std::vector<SMITH::IndexRange>> Get_Bagel_IndexRanges( std::shared_ptr<std::vector<std::string>> ranges_str);
      
      std::string get_det_name(std::shared_ptr<const Determinants> Detspace ) ;
       
      void get_civector_indexranges(int nstates) ;
      
      std::shared_ptr<SMITH::Tensor_<DataType>> convert_civec_to_tensor( std::shared_ptr<const Civec> civector, int state_num ) const;

      void build_detspace(std::shared_ptr<CIVecInfo<DataType>>  ci_info );

      void tester();

      
#ifndef NDEBUG

      bool gamma_2idx_contract_test( std::string gamma_name ) ;

      bool gamma_4idx_contract_test( std::string gamma_name ) ;

      bool build_sigma_4idx_tensor_tests(std::shared_ptr<GammaInfo<DataType>> gamma_4idx_info ) ;

#endif

};

}
#endif
