#ifndef __SRC_PROP_PROPTOOL_TENS_AND_CI_LIB_GAMMA_COMPUTER_H
#define __SRC_PROP_PROPTOOL_TENS_AND_CI_LIB_GAMMA_COMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/ci/fci/determinants.h>
#include <src/util/f77.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/vector_bundle.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

namespace bagel {

namespace Gamma_Computer { 
template<typename DataType>
class Gamma_Computer { 
 
    private: 

      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<const Determinants>>> bagel_determinant_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;

      std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<DataType>> tensor_calc_;
 
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<Vector_Bundle<DataType>>>>  sigma_data_map_;

      double thresh_;
      size_t civec_maxtile_;

    public: 
      Gamma_Computer(); 
      ~Gamma_Computer(){};

      ////////////////////// Setting map routines ////////////////////////
 
      void set_range_conversion_map (  std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map ) {
           range_conversion_map_ = range_conversion_map; return; }
 
      void set_gamma_info_map ( std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>>  gamma_info_map  ) {
           gamma_info_map_ = gamma_info_map; return; }

      void set_civec_data_map ( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map ) {
           civec_data_map_ = civec_data_map; return; }

      /////////////Gamma routines (Vector Bundle based) //////////////////
      
      void get_gamma( std::string gamma_name );

      void get_gammaN_from_sigmaN( std::shared_ptr<GammaInfo_Base> gamma_n_info );

      void compute_gammaN( std::shared_ptr<GammaInfo_Base> gammaN_info );
     
      void compute_sigmaN( std::shared_ptr<GammaInfo_Base> gammaN_info );

      void sigma_aa( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );
   
      void sigma2_aa( std::shared_ptr<Vector_Bundle<DataType>> sigma_aa, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                      std::shared_ptr<const Determinants> bra_det, std::shared_ptr<const Determinants> ket_det );

      void compute_sigma2( std::shared_ptr<GammaInfo_Base> gamma2_info );

      void sigma_bb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ); 

      void sigma2_bb( std::shared_ptr<Vector_Bundle<DataType>> sigma_bb, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<const Determinants> bra_det, std::shared_ptr<const Determinants> ket_det );

      void sigma_ab( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma2_ab( std::shared_ptr<Vector_Bundle<DataType>> sigma_ba, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<const Determinants> bra_det, std::shared_ptr<const Determinants> ket_det );

      void sigma_ba( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma2_ba( std::shared_ptr<Vector_Bundle<DataType>> sigma_ba, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                         std::shared_ptr<const Determinants> bra_det, std::shared_ptr<const Determinants> ket_det );

      void compute_eiej_on_ket( std::shared_ptr<Vector_Bundle<DataType>> eiej_on_ket, std::shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                std::shared_ptr<const Determinants> bra_det, std::shared_ptr<const Determinants> ket_det, std::string transition_type );

      /////////// Utility routines /////////////////////////
      
      std::vector<SMITH::IndexRange> get_bagel_indexranges( const std::vector<std::string>& ranges_str);

  
      void fill_and_link_determinant_map( int nele , int norb  );

      /////////// Variable access /////////////////////////
      
      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>>& gamma_info_map() { return gamma_info_map_ ; };
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>& gamma_data_map() {return gamma_data_map_; }
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>>& civec_data_map() {return civec_data_map_; }

      std::shared_ptr<SMITH::Tensor_<DataType>> gamma_data( std::string name ) {return gamma_data_map_->at(name); }
      std::shared_ptr<SMITH::Tensor_<DataType>> civec_data( std::string name ) {return civec_data_map_->at(name); }
};

}
}
#endif
