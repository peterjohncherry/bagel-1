#include <bagel_config.h>
#include <src/prop/proptool/moint_init.h>

using namespace std;
using namespace bagel;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
MOInt_Init<double>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref_,
                            int ncore, int nfrozenvirt, bool block_diag_fock ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   geom_ = geom; 
   //ref_ = ref;
//   ref = dynamic_pointer_cast<const RelReference>(ref);  
  shared_ptr<const Reference> bob = ref_ ; 
  ncore_ = ncore;
  nfrozenvirt_ = nfrozenvirt;
  block_diag_fock = block_diag_fock;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
MOInt_Init<std::complex<double>>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const RelReference> ref_,
                                              int ncore, int nfrozenvirt, bool block_diag_fock ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   geom_ = geom; 
   //ref_ = ref;
//   ref = dynamic_pointer_cast<const RelReference>(ref);  
  shared_ptr<const RelReference> bob = ref_ ; 
  ncore_ = ncore;
  nfrozenvirt_ = nfrozenvirt;
  block_diag_fock = block_diag_fock;
}




shared_ptr<const RDM<1>> MOInt_Init<double>::rdm1_av() const {
  return ref_->rdm1_av();
}

shared_ptr<const ZRDM<1>> MOInt_Init<complex<double>>::rdm1_av() const {
  return nullptr;
}


 // Need to do rel version
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//template<>
//bool MOInt_Init<double>::breit() const {
//      throw std::logic_error("Checking if Breit enabled in non-rel; this should not happen" ); return false; };

//template<>
//bool MOInt_Init<double>::gaunt() const {
//      throw std::logic_error("Checking if Gaunt enabled in non-rel; this should not happen" ); return false; };

//template<>
//bool MOInt_Init<std::complex<double>>::breit() const {return gaunt_; }

//template<>
//bool MOInt_Init<std::complex<double>>::gaunt() const {return breit_; }

//template<>
shared_ptr<const Matrix> MOInt_Init<double>::coeff() const {  return ref_->coeff();}


////////////////////////////////////////////////////////////////
//template<>
shared_ptr<const ZMatrix> MOInt_Init<complex<double>>::coeff() const {
////////////////////////////////////////////////////////////////
 //shared_ptr<const ZCoeff_Striped> c ;//= dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
//  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
    shared_ptr<const ZMatrix> bob; return bob;
}
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////
template class MOInt_Init<double>;
template class MOInt_Init<std::complex<double>>;
////////////////////////////////////////////////////////////////
//                                 geom_(geom), ref_(ref), ncore_(ncore), nfrozenvirt_(nfrozenvirt),
//
//                                 block_diag_fock_(block_diag_fock) {}
