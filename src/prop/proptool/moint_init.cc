
#include <bagel_config.h>
#include <src/prop/proptool/moint_init.h>

using namespace std;
using namespace bagel;

template<typename DataType> 
MOInt_Init<DataType>::MOInt_Init(std::shared_ptr<const Geometry> geom,  std::shared_ptr<const Reference> ref_ ) :
                                 geom_(geom), ref_(ref) {} 

template<>
shared_ptr<const Matrix> MOInt_Init<double>::coeff() const {
  return ref_->coeff();
}

template<>
shared_ptr<const ZMatrix> MOInt_Init<complex<double>>::coeff() const {
  shared_ptr<const ZCoeff_Striped> c = dynamic_pointer_cast<const RelReference>(ref_)->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
}

