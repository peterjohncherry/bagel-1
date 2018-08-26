#ifndef __SRC_PROPTOOL_MOINT_COMPUTER_H
#define __SRC_PROPTOOL_MOINT_COMPUTER_H


// Should presumably be expanded. Written like this to keep things seperated. 
// Should be done so all calls to members of info_ from K2ext_new and MOFock_new just refer
// to members of MO_integrator defined from info_. So K2ext_new and MOFock_new might work best
// as sub classes of MO_Integrator, but I think this might be a headache with the templates, 
// so I will leave it for now.

#include <stddef.h>
#include <memory>
#include <stdexcept>
#include <complex>
#include <src/prop/proptool/integrals/moint_init.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>

namespace bagel {
  template<typename DataType>
  class MOInt_Computer {
    private :
      using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

      using Tensor = typename std::conditional<std::is_same<DataType,double>::value,
                                               SMITH::Tensor_<double>,SMITH::Tensor_<std::complex<double>>>::type;

      std::shared_ptr<const MOInt_Init<DataType>> info_;
      std::shared_ptr<const MatType> coeffs_;
      std::shared_ptr<std::map<std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_;

      std::shared_ptr<SMITH::Tensor_<DataType>> f1_;
      std::shared_ptr<SMITH::Tensor_<DataType>> h1_;
      std::shared_ptr<SMITH::Tensor_<DataType>> v2_;
      
      bool got_fock_coeffs_; 

    public :
      MOInt_Computer( std::shared_ptr<const MOInt_Init<DataType>> r,
                      std::shared_ptr<std::map<std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map  ):
                      info_(r), coeffs_(info_->coeff()), range_conversion_map_(range_conversion_map)  {}
      ~MOInt_Computer(){};

      //note, this does not have the diagonal component
      void calculate_v2( const std::vector<SMITH::IndexRange>& blocks );
      void calculate_v2( const std::vector<std::string>& blocks );

      //is the core fock minus diagonal component from above
      void calculate_h1( const std::vector<std::string>& blocks, bool set_coeffs = false, bool set_fock = false );
      void calculate_h1( const std::vector<SMITH::IndexRange>& blocks, bool set_coeffs = false, bool set_fock = false );
 
      // core fock, without subtracted diagonal component
      void calculate_fock( const std::vector<std::string>& blocks, bool set_coeffs = false, bool set_h1 = true );
      void calculate_fock( const std::vector<SMITH::IndexRange>& blocks, bool set_coeffs = false, bool set_h1 = true );
     
      std::shared_ptr<const MatType> coeffs() { return coeffs_; } 
      std::shared_ptr<SMITH::Tensor_<DataType>> f1() { return f1_; } 
      std::shared_ptr<SMITH::Tensor_<DataType>> h1() { return h1_; } 
      std::shared_ptr<SMITH::Tensor_<DataType>> v2() { return v2_; } 

      bool got_fock_coeffs() const { return got_fock_coeffs_; }

      //TEST
      std::shared_ptr<SMITH::Tensor_<DataType>> t_from_smith_;
      std::shared_ptr<SMITH::Tensor_<DataType>> s_from_smith_;
      std::shared_ptr<SMITH::Tensor_<DataType>> v2_from_smith_;
      std::shared_ptr<SMITH::Tensor_<DataType>> v2_smith_;
      std::shared_ptr<SMITH::Tensor_<DataType>> build_s_test_tensor( const std::vector<int>& ordering);
      std::shared_ptr<SMITH::Tensor_<DataType>> calculate_v2_smith();
      std::shared_ptr<SMITH::Tensor_<DataType>> get_test_tensor( const std::vector<SMITH::IndexRange>& blocks  );
      std::shared_ptr<SMITH::Tensor_<DataType>> get_test_tensor( const std::vector<std::string>& blocks  );
      // ENDTEST

  };
}
#endif
