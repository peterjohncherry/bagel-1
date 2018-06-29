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

namespace B_Gamma_Computer { 
template<typename DataType>
class B_Gamma_Computer { 
 
    private: 

      std::shared_ptr<const Dvec> cc_; 
      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map_;

      std::shared_ptr<std::map< std::string, std::shared_ptr<Dvec>>> dvec_sigma_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<Civec>>> cvec_old_map_; 
      std::shared_ptr<std::map< std::string, std::shared_ptr<Determinants>>> det_old_map_;

      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;
      
      std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<DataType>> tensor_calc_;
 
      std::shared_ptr<std::map< std::string, std::shared_ptr<Vector_Bundle<DataType>>>> new_sigma_data_map_;

    public: 

      B_Gamma_Computer( std::shared_ptr<const Dvec> cc_in ); 
      ~B_Gamma_Computer(){};
      
      void set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map__in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map__in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map,  
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map );
     
      ////////////Gamma routines (RDM class based) //////////////////
      
      void get_gamma( std::string gamma_name );
      
      void get_wfn_data( std::shared_ptr<CIVecInfo_Base>  cvec_info );
      
      void compute_sigma2( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void compute_sigmaN( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void get_gamma2_from_sigma2( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void get_gammaN_from_sigmaN( std::shared_ptr<GammaInfo_Base> gammaN_info );
      
      void sigma_2a1( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<Determinants> dets );
      
      void sigma_2a2( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<Determinants> dets ) ;
      
      void sigma_2ab( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<Determinants> bra_det, std::shared_ptr<Determinants> ket_det );

      void fill_and_link_determinant_map( int nelec, int norb ); 

      void sigma_ab( std::shared_ptr<GammaInfo_Base> gamma_info ) ; 

      void sigma_aa_test( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma = true );

      void sigma_bb_test( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma = false );

      void sigma2_test( std::shared_ptr<GammaInfo_Base> gamma_info );

      void sigma_ab_test( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma_ba_test( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );

      void sigma_aa_vb( std::shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma );
   
      void compute_sigma2_vb( std::shared_ptr<GammaInfo_Base> gamma2_info );

      /////////// Utility routines /////////////////////////
      
      std::shared_ptr<std::vector<SMITH::IndexRange>>  Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);
      
      void convert_Dvec_sigma_to_tensor( std::shared_ptr<GammaInfo_Base> gammaN_info );
      
      void convert_civec_to_tensor( std::string civec_name ) ;
      
      /////////// Varaible access /////////////////////////
      
      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map() { return gamma_info_map_ ; };
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map() {return gamma_data_map_; }
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map() {return sigma_data_map_; }
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map() {return civec_data_map_; }

      std::shared_ptr<SMITH::Tensor_<DataType>> gamma_data( std::string name ) {return gamma_data_map_->at(name); }
      std::shared_ptr<SMITH::Tensor_<DataType>> sigma_data( std::string name ) {return sigma_data_map_->at(name); }
      std::shared_ptr<SMITH::Tensor_<DataType>> civec_data( std::string name ) {return civec_data_map_->at(name); }
};

}
}
#endif
