
#include <bagel_config.h>
#include <src/prop/proptool/moint_init.h>

using namespace std;
using namespace bagel;

template<typename DataType> 
MOInt_Init<DataType>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref_,
                                 int ncore, int nfrozenvirt, bool block_diag_fock ):            
                                 geom_(geom), ref_(ref), ncore_(ncore), nfrozenvirt_(nfrozenvirt),
                                 block_diag_fock_(block_diag_fock) {}
 // Need to do rel version

template<>
bool MOInt_Init<double>::breit() const {
      throw std::logic_error("Checking if Breit enabled in non-rel; this should not happen" ); return false; };

template<>
bool MOInt_Init<double>::gaunt() const {
      throw std::logic_error("Checking if Gaunt enabled in non-rel; this should not happen" ); return false; };

template<>
bool MOInt_Init<std::complex<double>>::breit() const {return gaunt_; }

template<>
bool MOInt_Init<std::complex<double>>::gaunt() const {return breit_; }

template<>
shared_ptr<const Matrix> MOInt_Init<double>::coeff() const {
  return ref_->coeff();
}

template<>
shared_ptr<const ZMatrix> MOInt_Init<complex<double>>::coeff() const {
  shared_ptr<const ZCoeff_Striped> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
}

