#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_B_GAMMA_COMPUTER_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_B_GAMMA_COMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/ci/fci/civec.h>
#include <src/ci/fci/dvec.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/vector_bundle.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

namespace bagel {

namespace Gamma_Computer_Init { 
template<typename DataType>
class Gamma_Computer_Init { 
  public: 

    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_; 

    Gamma_Computer_Init( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map ): 
                         range_conversion_map_(range_conversion_map),
                         civec_data_map_( std::make_shared<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>()) {};
    ~Gamma_Computer_Init(){};
    
     void add_civec( std::shared_ptr<Civec> civec, int state_num );
    
   
 };
};
namespace Gamma_Computer { 
template<typename DataType>
class Gamma_Computer { 
 
    private: 

      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map_;

      std::shared_ptr<std::map< std::string, std::shared_ptr<Determinants>>> bagel_determinant_map_;

      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;
      
      std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<DataType>> tensor_calc_;
 
      std::shared_ptr<std::map< std::string, std::shared_ptr<Vector_Bundle<DataType>>>>  vb_sigma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> vb_gamma_data_map_;

      double thresh_;
      size_t civec_maxtile_;

    public: 

      Gamma_Computer(); 
      Gamma_Computer( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map,
                      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map );
      ~Gamma_Computer(){};
      

      ////////////Gamma routines (Vector Bundle based) //////////////////
      
      void get_gamma_vb( std::string gamma_name );

      void get_gammaN_from_sigmaN_vb( std::shared_ptr<GammaInfo_Base> gamma_n_info );

      void compute_gammaN_vb( std::shared_ptr<GammaInfo_Base> gammaN_info );
     
      void compute_sigmaN_vb( std::shared_ptr<GammaInfo_Base> gammaN_info );

      void sigma_aa_vb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );
   
      void sigma2_aa_vb( std::shared_ptr<Vector_Bundle<DataType>> sigma_aa, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det );

      void compute_sigma2_vb( std::shared_ptr<GammaInfo_Base> gamma2_info );

      void sigma_bb_vb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ); 

      void sigma2_bb_vb( std::shared_ptr<Vector_Bundle<DataType>> sigma_bb, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det );

      void sigma_ab_vb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma2_ab_vb( std::shared_ptr<Vector_Bundle<DataType>> sigma_ba, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det );

      void sigma_ba_vb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma2_ba_vb( std::shared_ptr<Vector_Bundle<DataType>> sigma_ba, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det );

      void compute_eiej_on_ket( std::shared_ptr<Vector_Bundle<DataType>> eiej_on_ket, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det, std::string transition_type );

      /////////// Utility routines /////////////////////////
      
      std::vector<SMITH::IndexRange>  get_bagel_indexranges( const std::vector<std::string>& ranges_str);
      void get_wfn_data( std::shared_ptr<CIVecInfo_Base>  cvec_info );
      void convert_civec_to_tensor( std::string civec_name ) ;
};

}
}
#endif
