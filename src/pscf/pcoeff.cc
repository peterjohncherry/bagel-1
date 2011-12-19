//
// Author : Toru Shiozaki
// Date   : July 2009
//

#include <src/pscf/pcoeff.h>
#include <src/pscf/f77.h>

typedef std::complex<double> Complex;

using namespace std;

PCoeff::PCoeff(const PMatrix1e& inp) : PMatrix1e(inp.geom(), inp.ndim(), inp.mdim()) {

  const int unit = 1;
  const int ndim_ = inp.ndim();
  const int mdim_ = inp.mdim();
  zcopy_(&totalsize_, inp.data()->front(), &unit, data()->front(), &unit); 

}


PCoeff::PCoeff(const shared_ptr<PGeometry> gm, const int ndim, const int mdim)
 : PMatrix1e(gm, ndim, mdim) {


}

PCoeff::~PCoeff() {

}

PMatrix1e PCoeff::form_density_rhf(const bool return_ao) const {
  const int nocc = geom_->nocc() / 2;
  assert(geom_->nocc() % 2 == 0);
  const Complex one(1.0, 0.0);
  const Complex zero(0.0, 0.0);

  // first, form the density matrix in k space. 
  PMatrix1e k_density(geom_);
  int kcount = 0;
  for (int k = -K(); k <= K(); ++k, ++kcount) { 
    const int koffset = kcount * blocksize_;
    zgemm_("N", "C", &nbasis_, &nbasis_, &nocc, &one, data()->pointer(koffset), &nbasis_, 
                                                      data()->pointer(koffset), &nbasis_,
                                     &zero, k_density.data()->pointer(koffset), &nbasis_); 
  }

  // back Fourier transform
  if (return_ao)
      return k_density.bft();
  else
      return k_density;
}


pair<shared_ptr<PCoeff>, shared_ptr<PCoeff> > PCoeff::split(const int nrow1, const int nrow2) {
  shared_ptr<PCoeff> out1(new PCoeff(geom_, nrow1, mdim_));
  shared_ptr<PCoeff> out2(new PCoeff(geom_, nrow2, mdim_));

  assert(nrow1+nrow2 == ndim_);
  assert(blocksize_ == out1->blocksize() + out2->blocksize());
  assert(blocksize_ == ndim_ * mdim_);

  Complex* source = data_->front();
  Complex* data1 = out1->data()->front();
  Complex* data2 = out2->data()->front();

  for (int i = -K(); i <= K(); ++i) {
    for (int m = 0; m != mdim_; ++m, data1+=out1->ndim(), data2+=out2->ndim(), source+=ndim_) {
      copy(source,       source+nrow1,       data1);
      copy(source+nrow1, source+nrow1+nrow2, data2);
    }
  }

  return make_pair(out1, out2);
}
