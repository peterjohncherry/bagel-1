#include <bagel_config.h>
#include <src/prop/proptool/moint_init.h>

using namespace std;
using namespace bagel;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
MOInt_Init<double>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref,
                            int ncore, int nfrozenvirt, bool block_diag_fock ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
   geom_ = geom; 
   ref_ = ref;
//   ref = dynamic_pointer_cast<const RelReference>(ref);  
  shared_ptr<const Reference> bob = ref_ ; 
  ncore_ = ncore;
  nfrozenvirt_ = nfrozenvirt;
  block_diag_fock = block_diag_fock;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
MOInt_Init<std::complex<double>>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const RelReference> ref,
                                              int ncore, int nfrozenvirt, bool block_diag_fock ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  geom_ = geom; 
  ref_ = ref;
  ref = dynamic_pointer_cast<const RelReference>(ref);  
  shared_ptr<const RelReference> bob = ref_ ; 
  ncore_ = ncore;
  nfrozenvirt_ = nfrozenvirt;
  block_diag_fock = block_diag_fock;
}

////////////////////////////////////////////////////////////////
//template<>
shared_ptr<const ZMatrix> MOInt_Init<complex<double>>::coeff() const {
////////////////////////////////////////////////////////////////
shared_ptr<const ZCoeff_Striped> c = ref_->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
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
