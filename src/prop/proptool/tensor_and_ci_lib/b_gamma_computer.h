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
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> new_gamma_data_map_;

      double thresh_;
      size_t civec_maxtile_;

    public: 

      B_Gamma_Computer( std::shared_ptr<const Dvec> cc_in ); 
      ~B_Gamma_Computer(){};
      
      void set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map__in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map__in,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> gamma_data_map,
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> sigma_data_map,  
                     std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<DataType>>>> civec_data_map );
     

      ////////////Gamma routines (Dvec based) //////////////////
      
      void get_gamma( std::string gamma_name );
      
      void compute_sigma2( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void compute_sigmaN( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void get_gamma2_from_sigma2( std::shared_ptr<GammaInfo_Base> gamma_info );
      
      void get_gammaN_from_sigmaN( std::shared_ptr<GammaInfo_Base> gammaN_info );
    
      void sigma_2a1( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<Determinants> dets );
      
      void sigma_2a2( DataType* cvec_ptr, DataType* sigma_ptr, std::shared_ptr<Determinants> dets ) ;

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
      
      std::shared_ptr<std::vector<SMITH::IndexRange>>  Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);
      
      void convert_Dvec_sigma_to_tensor( std::shared_ptr<GammaInfo_Base> gammaN_info );
      
      void get_wfn_data( std::shared_ptr<CIVecInfo_Base>  cvec_info );

      void convert_civec_to_tensor( std::string civec_name ) ;
  
      void fill_and_link_determinant_map( int nele , int norb  );

      /////////// Variable access /////////////////////////
      
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
