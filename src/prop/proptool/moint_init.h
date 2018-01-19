#ifndef __SRC_PROPTOOL_MOINT_INFO_H
#define __SRC_PROPTOOL_MOINT_INFO_H

#include <stddef.h>
#include <memory>
#include <stdexcept>
#include <src/wfn/relreference.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>

namespace bagel {

  template<typename DataType>
  class MOInt_Init {
    private :
      using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

      using Tensor = typename std::conditional<std::is_same<DataType,double>::value,
                                               SMITH::Tensor_<double>,SMITH::Tensor_<std::complex<double>>>::type;

      using IndexRange = SMITH::IndexRange;

      int nclosed_;
      int nocc_;
      int nact_;
      int ncore_;
      int nvirt_;
      int nfrozenvirt_;

      bool block_diag_fock_;
      bool breit_;
      bool gaunt_;

      std::string method_;

      std::shared_ptr<const Geometry> geom_;
      std::shared_ptr<const Reference> ref_;

    public:
      MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref_ );
      ~MOInt_Init(){};

      int nclosed()const{ return nclosed_;}
      int nocc()const { return nocc_;}
      int nact()const{ return nact_;}
      int ncore()const{ return ncore_;}
      int nvirt()const { return nvirt_;}
      int nfrozenvirt()const{ return nfrozenvirt_;}
      
      bool block_diag_fock()const { return block_diag_fock_; }     
      bool breit()const { return breit_; }     
      bool gaunt()const { return gaunt_; }     

      std::string method()const{ return method_; }

      std::shared_ptr<const Geometry> geom() const { return geom_;}
      std::shared_ptr<const RDM<1>> rdm1_av() const { return  ref_->rdm1_av();}
      std::shared_ptr<const Hcore> hcore() const { return ref_->hcore();}
      std::shared_ptr<const MatType> coeff() const { assert(false);}
      std::shared_ptr<const Reference> ref() const { return ref_ ; }

  };

template<> std::shared_ptr<const Matrix> MOInt_Init<double>::coeff() const;
template<> std::shared_ptr<const ZMatrix> MOInt_Init<std::complex<double>>::coeff() const;

//extern template class MOInt_Init<double>;
//extern template class MOInt_Init<std::complex<double>>;

};
#endif
