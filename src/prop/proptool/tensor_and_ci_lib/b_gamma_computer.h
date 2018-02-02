#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_B_GAMMA_COMPUTER_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_B_GAMMA_COMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/ci/fci/civec.h>
#include <src/ci/fci/dvec.h>
#include <src/util/f77.h>
#include <src/prop/proptool/algebraic_manipulator/expression.h> 
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>

namespace bagel {

namespace B_Gamma_Computer { 

class B_Gamma_Computer { 

    public: 
    B_Gamma_Computer( std::shared_ptr<const Dvec> cc_in, 
                       std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_in,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>>  Gamma_info_map_in ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Gamma_data_map_in ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Sigma_data_map_in ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<SMITH::Tensor_<double>>>> CIvec_data_map_in );
    ~B_Gamma_Computer(){};
  
    std::shared_ptr<const Dvec> cc_; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> Gamma_info_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Gamma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> Sigma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> CIvec_data_map;

    std::shared_ptr<std::map< std::string, std::shared_ptr<Dvec>>> dvec_sigma_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Civec>>> cvec_old_map; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<Determinants>>> det_old_map;

    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map;

    std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> Tensor_Calc;

    ////////////Gamma routines (RDM class based) //////////////////

    void get_gamma( std::string gamma_name );
    
    void get_wfn_data( std::shared_ptr<CIVecInfo<double>>  cvec_info );

    void compute_sigma2( std::shared_ptr<GammaInfo> gamma_info );

    void compute_sigmaN( std::shared_ptr<GammaInfo> gamma_info );
  
    void get_gamma2_from_sigma2( std::shared_ptr<GammaInfo> gamma_info );

    void get_gammaN_from_sigmaN( std::shared_ptr<GammaInfo> gammaN_info );

    void sigma_2a1( double* cvec_ptr, double* sigma_ptr, std::shared_ptr<Determinants> dets );

    void sigma_2a2( double* cvec_ptr, double* sigma_ptr, std::shared_ptr<Determinants> dets ) ;

    /////////// Utility routines /////////////////////////

    std::shared_ptr<std::vector<SMITH::IndexRange>>  Get_Bagel_IndexRanges(std::shared_ptr<std::vector<std::string>> ranges_str);

    void convert_Dvec_sigma_to_tensor( std::shared_ptr<GammaInfo> gammaN_info );

    void convert_civec_to_tensor( std::string civec_name ) ;

};

}
}
#endif
