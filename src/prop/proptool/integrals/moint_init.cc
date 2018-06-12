#include <bagel_config.h>
#include <src/prop/proptool/integrals/moint_init.h>

using namespace std;
using namespace bagel;

////////////////////////////////////////////////////////////////////////
//TODO Replace this with something proper relativistic version
shared_ptr<const ZMatrix> MOInt_Init<complex<double>>::coeff() const {
////////////////////////////////////////////////////////////////////////
shared_ptr<const ZCoeff_Striped> c = ref_->relcoeff();
  return c->block_format(nclosed(), nact(), nvirt()+nfrozenvirt(), 0);
    shared_ptr<const ZMatrix> bob; return bob;
}
/////////////////////////////////////////////////////////////////////////
template class MOInt_Init<double>;
template class MOInt_Init<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
