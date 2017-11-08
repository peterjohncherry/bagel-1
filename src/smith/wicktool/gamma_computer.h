#ifndef __SRC_SMITH_WICKTOOL_GAMMA_COMPUTER_H
#define __SRC_SMITH_WICKTOOL_GAMMA_COMPUTER_H

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/wicktool/equation_computer.h>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic.h>
#include <src/smith/wicktool/tensor_arithmetic_utils.h>
#include <src/smith/wicktool/states_info.h>

namespace bagel {
namespace SMITH { 

namespace Gamma_Computer { 

class Gamma_Computer { 

    public: 
      Gamma_Computer( std::shared_ptr< std::map< std::string, std::shared_ptr<GammaInfo>>>          GammaMap,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<Tensor_<double>>>>    CIvec_data_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<Tensor_<double>>>>    Sigma_data_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<Tensor_<double>>>>    Gamma_data_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<const Determinants>>> Determinants_map,
                      std::shared_ptr< std::map< std::string, std::shared_ptr<IndexRange>>>         range_conversion_map ); 
      ~Gamma_Computer(){};
      
      int nelea_  ;
      int neleb_  ;
      int norb_   ;
      int nstate_ ;
      int maxtile ;
      int cimaxblock;
      std::shared_ptr<IndexRange> virt_   ;
      std::shared_ptr<IndexRange> active_ ;
      std::shared_ptr<IndexRange> closed_ ;
      
      std::shared_ptr<const Dvec> cc_; 
      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Gamma_data_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Sigma_data_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CIvec_data_map;

      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> GammaMap;
      std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;
      std::shared_ptr<std::map< std::string, std::shared_ptr<const Determinants>>> Determinants_map;
      
      std::shared_ptr<Equation<double>> Eqn_info;
      
      std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> Tensor_Calc;
      
      /////////// Gamma routines /////////////////////////
      void get_gamma_tensor( std::string gamma_name);
      
      void build_sigma_2idx_tensor( std::shared_ptr<GammaInfo> gamma_info ) ;
      
      void build_gamma_2idx_tensor( std::string gamma_name ) ;

      void build_gamma_4idx_tensor( std::string gamma_name ) ;
      
      void build_sigma_block( std::shared_ptr<Tensor_<double>> sigma_tens, std::vector<Index>& sigma_id_blocks, std::vector<int>& sigma_offsets, std::string Ket_name  ) const ;
      
      void get_civector_indexranges(int nstates) ;
      
      std::shared_ptr<Tensor_<double>> convert_civec_to_tensor( std::shared_ptr<const Civec> civector, int state_num ) const;
      
      /////////// Utility routines /////////////////////////
      
      std::shared_ptr<std::vector<IndexRange>> Get_Bagel_IndexRanges( std::shared_ptr<std::vector<std::string>> ranges_str);
      
      std::string get_sigma_name( std::string Bra_name, std::string Ket_name , std::shared_ptr<std::vector<std::string>>  orb_ranges, std::shared_ptr<std::vector<bool>> aops ) ;
      
      std::string get_det_name(std::shared_ptr<const Determinants> Detspace ) ;
      
      void tester();

};

}
}
}
#endif
