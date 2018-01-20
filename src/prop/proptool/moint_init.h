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

      int ncore_;
      int nfrozenvirt_;

      bool block_diag_fock_;
      bool breit_;
      bool gaunt_;

      std::string method_;

      std::shared_ptr<const Geometry> geom_;
      std::shared_ptr<const Reference> ref_;

    public:
      MOInt_Init( std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref, 
                  int ncore, int nfrozenvirt, bool block_diag_fock );

      ~MOInt_Init(){};

      int ncore() const { return ncore_; }
      int nclosed() const { return ref_->nclosed(); }
      int nact() const { return ref_->nact(); }
      int nocc() const { return ref_->nocc(); }
      int nvirt() const { return ref_->nvirt() - nfrozenvirt_; }
      int nfrozenvirt() const { return nfrozenvirt_; }
      
      bool block_diag_fock()const { return block_diag_fock_; }     
      bool breit()const ;
      bool gaunt()const ;

      std::string method()const{ return method_; }

      std::shared_ptr<const Geometry> geom() const { return geom_;}
      std::shared_ptr<const RDM<1>> rdm1_av() const { return  ref_->rdm1_av();}
      std::shared_ptr<const Hcore> hcore() const { return ref_->hcore();}
      std::shared_ptr<const MatType> coeff() const  {  return ref_->coeff();}
      std::shared_ptr<const Reference> ref() const { return ref_ ; }

  };

  template<> class MOInt_Init<double> {
    private :

      std::shared_ptr<const Geometry> geom_;
      std::shared_ptr<const Reference> ref_;

      int ncore_;
      int nfrozenvirt_;

      bool block_diag_fock_;
      std::string method_;

    public:
      MOInt_Init( std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref, 
                  int ncore, int nfrozenvirt, bool block_diag_fock );

      ~MOInt_Init(){};

      int ncore() const { return ncore_; }
      int nclosed() const { return ref_->nclosed(); }
      int nact() const { return ref_->nact(); }
      int nocc() const { return ref_->nocc(); }
      int nvirt() const { return ref_->nvirt() - nfrozenvirt_; }
      int nfrozenvirt() const { return nfrozenvirt_; }
      
      bool block_diag_fock()const { return block_diag_fock_; }     
      bool breit() const { return false; }     
      bool gaunt() const { return false; }     

      std::string method()const{ return method_; }

      std::shared_ptr<const Geometry> geom() const { return geom_;}
      std::shared_ptr<const RDM<1>> rdm1_av() const { return ref_->rdm1_av(); } 
      std::shared_ptr<const Hcore> hcore() const { return ref_->hcore();}
      std::shared_ptr<const Matrix> coeff() const { return ref_->coeff(); }
      std::shared_ptr<const Reference> ref() const { return ref_ ; }

  };

  template<> class MOInt_Init<std::complex<double>> {
    private :

      std::shared_ptr<const Geometry> geom_;
      std::shared_ptr<const RelReference> ref_;

      int ncore_;
      int nfrozenvirt_;

      bool block_diag_fock_;
      bool breit_;
      bool gaunt_;

      std::string method_;

    public:
      MOInt_Init( std::shared_ptr<const Geometry> geom,  std::shared_ptr<const RelReference> ref, 
                  int ncore, int nfrozenvirt, bool block_diag_fock );

      ~MOInt_Init(){};

      int ncore() const { return ncore_; }
      int nclosed() const { return ref_->nclosed(); }
      int nact() const { return ref_->nact(); }
      int nocc() const { return ref_->nocc(); }
      int nvirt() const { return ref_->nvirt() - nfrozenvirt_; }
      int nfrozenvirt() const { return nfrozenvirt_; }
      
      bool block_diag_fock()const { return block_diag_fock_; }     
      bool breit()const { ref_->breit();}     
      bool gaunt()const { ref_->gaunt();}     

      std::string method()const{ return method_; }

      std::shared_ptr<const ZRDM<1>> rdm1_av() const { return nullptr; } ; // weird, but this is how it is defined in smith_info
      std::shared_ptr<const ZMatrix> coeff() const;

      std::shared_ptr<const Geometry> geom() const { return geom_;}
      std::shared_ptr<const Hcore> hcore() const { return ref_->hcore();}
      std::shared_ptr<const RelReference> ref() const { return ref_ ; }

  };


//template<> std::shared_ptr<const Matrix> MOInt_Init<double>::coeff() const;
//template<> std::shared_ptr<const ZMatrix> MOInt_Init<std::complex<double>>::coeff() const;

//template<> bool MOInt_Init<double>::breit() const;
//template<> bool MOInt_Init<double>::gaunt() const;
//template<> bool MOInt_Init<std::complex<double>>::breit()const ;
//template<> bool MOInt_Init<std::complex<double>>::gaunt()const ;

};
#endif
