//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: denom.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/denom.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace SMITH;
using namespace btas;

template<typename DataType>
Denom<DataType>::Denom(shared_ptr<const MatType> fock, const int nstates, const double th) : fock_(fock), thresh_(th) {
  const size_t ndim = fock->mdim() * nstates;
  const size_t ndim2 = fock->mdim() * ndim;
  const size_t ndim3 = fock->mdim() * ndim2;
  const int fac2 = is_same<DataType,double>::value ? 2 : 1;

  shalf_x_ = make_shared<MatType>(ndim, ndim);
  shalf_h_ = make_shared<MatType>(ndim, ndim);
  shalf_xx_ = make_shared<MatType>(ndim2, ndim2);
  shalf_hh_ = make_shared<MatType>(ndim2, ndim2);
  shalf_xh_ = make_shared<MatType>(ndim2*fac2, ndim2*fac2);
  shalf_xxh_ = make_shared<MatType>(ndim3, ndim3);
  shalf_xhh_ = make_shared<MatType>(ndim3, ndim3);

  work_x_ = make_shared<MatType>(ndim, ndim);
  work_h_ = make_shared<MatType>(ndim, ndim);
  work_xx_ = make_shared<MatType>(ndim2, ndim2);
  work_hh_ = make_shared<MatType>(ndim2, ndim2);
  work_xh_ = make_shared<MatType>(ndim2*fac2, ndim2*fac2);
  work_xxh_ = make_shared<MatType>(ndim3, ndim3);
  work_xhh_ = make_shared<MatType>(ndim3, ndim3);

  denom_x_ = VectorB(ndim);
  denom_h_ = VectorB(ndim);
  denom_xx_ = VectorB(ndim2);
  denom_hh_ = VectorB(ndim2);
  denom_xh_ = VectorB(ndim2*fac2);
  denom_xxh_ = VectorB(ndim3);
  denom_xhh_ = VectorB(ndim3);
}


template<typename DataType>
void Denom<DataType>::append(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                           shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  init_x_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_h_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_xx_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_hh_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_xh_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_xxh_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
  init_xhh_(jst, ist, rdm1, rdm2, rdm3, rdm4f);
}


template<typename DataType>
void Denom<DataType>::compute() {
  // TODO in principle, we can use smaller dimension (i.e., use canonical orthogonalization for shalf_x_)
  {
    shalf_x_->inverse_half(thresh_);
    MatType tmp(*shalf_x_ % *work_x_ * *shalf_x_);
    tmp.diagonalize(denom_x_);
    shalf_x_ = make_shared<MatType>(tmp % *shalf_x_);
  }
  {
    shalf_h_->inverse_half(thresh_);
    MatType tmp(*shalf_h_ % *work_h_ * *shalf_h_);
    tmp.diagonalize(denom_h_);
    shalf_h_ = make_shared<MatType>(tmp % *shalf_h_);
  }
  {
    shalf_xx_->inverse_half(thresh_);
    MatType tmp(*shalf_xx_ % *work_xx_ * *shalf_xx_);
    tmp.diagonalize(denom_xx_);
    shalf_xx_ = make_shared<MatType>(tmp % *shalf_xx_);
  }
  {
    shalf_hh_->inverse_half(thresh_);
    MatType tmp(*shalf_hh_ % *work_hh_ * *shalf_hh_);
    tmp.diagonalize(denom_hh_);
    shalf_hh_ = make_shared<MatType>(tmp % *shalf_hh_);
  }
  {
    shalf_xh_->inverse_half(thresh_);
    MatType tmp(*shalf_xh_ % *work_xh_ * *shalf_xh_);
    tmp.diagonalize(denom_xh_);
    shalf_xh_ = make_shared<MatType>(tmp % *shalf_xh_);
  }
  {
    shalf_xxh_->inverse_half(thresh_);
    MatType tmp(*shalf_xxh_ % *work_xxh_ * *shalf_xxh_);
    tmp.diagonalize(denom_xxh_);
    shalf_xxh_ = make_shared<MatType>(tmp % *shalf_xxh_);
  }
  {
    shalf_xhh_->inverse_half(thresh_);
    MatType tmp(*shalf_xhh_ % *work_xhh_ * *shalf_xhh_);
    tmp.diagonalize(denom_xhh_);
    shalf_xhh_ = make_shared<MatType>(tmp % *shalf_xhh_);
  }
}


template<typename DataType>
void Denom<DataType>::init_x_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                            shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim = nact;

  shalf_x_->copy_block(dim*jst, dim*ist, dim, dim, rdm1);

  MatType work2(dim, dim);
  auto rdm2v = group(group(*rdm2, 2,4),0,2);
  auto workv = group(work2, 0,2);
  contract(1.0, rdm2v, {0,1}, group(*fock_,0,2), {1}, 0.0, workv, {0});
  work_x_->copy_block(dim*jst, dim*ist, dim, dim, work2);
}


template<typename DataType>
void Denom<DataType>::init_h_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                            shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact;
  auto shalf = make_shared<MatType>(dim, dim);
  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  if (jst == ist) {
    copy_n(rdm1->data(), rdm1->size(), shalf->data());
    if (is_same<DataType,complex<double>>::value)
      shalf = shalf->get_conjg();
    shalf->scale(-1.0);
    shalf->add_diag(fac2); //.. making hole 1RDM
  } else {
    blas::transpose(rdm1->data(), dim, dim, shalf->data());
    shalf->scale(-1.0);
  }
  shalf_h_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<2,DataType>> ovl = rdm2->clone();
  sort_indices<1,0,2,0,1,-1,1>(rdm2->data(), ovl->data(), nact, nact, nact*nact);

  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        // see Celani eq. A7
        ovl->element(k, i, j, i) +=        shalf->element(k,j);
        ovl->element(i, k, i, j) += -1.0  * rdm1->element(k,j);
        ovl->element(i, i, k, j) +=  fac2 * rdm1->element(k,j);
      }
  MatType work2(dim, dim);
  auto ovlv = group(group(*ovl,2,4),0,2);
  auto workv = group(work2, 0,2);
  contract(1.0, ovlv, {0,1}, group(*fock_,0,2), {1}, 0.0, workv, {0});
  work_h_->copy_block(dim*jst, dim*ist, dim, dim, work2);
}



template<typename DataType>
void Denom<DataType>::init_xx_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                             shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  auto shalf = make_shared<MatType>(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(rdm2->data(), shalf->data(), nact, nact, nact, nact);
  shalf_xx_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  auto work2 = make_shared<MatType>(dim, dim);
  auto workv = group(*work2, 0,2);
  auto rdm3v = group(group(*rdm3,4,6),0,4);
  contract(1.0, rdm3v, {0,1},  group(*fock_,0,2), {1}, 0.0, workv, {0});
  auto work = make_shared<MatType>(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work->data(), nact, nact, nact, nact);
  work_xx_->copy_block(dim*jst, dim*ist, dim, dim, work);
}


template<typename DataType>
void Denom<DataType>::init_hh_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                             shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  shared_ptr<RDM<2,DataType>> shalf = rdm2->clone();
  sort_indices<1,3,0,2,0,1,1,1>(rdm2->data(), shalf->data(), nact, nact, nact, nact);
  const double fac = jst == ist ? 1.0 : 0.0;
  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  const double fac4 = is_same<DataType,double>::value ? 4.0 : 1.0;
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i1 = 0; i1 != nact; ++i1) {
      shalf->element(i2, i1, i2, i1) +=  fac4 * fac;
      shalf->element(i2, i1, i1, i2) += -fac2 * fac;
      for (int i0 = 0; i0 != nact; ++i0) {
        shalf->element(i0, i1, i1, i2) +=         rdm1->element(i2, i0);
        shalf->element(i0, i1, i2, i1) += -fac2 * rdm1->element(i2, i0);
        shalf->element(i1, i0, i1, i2) += -fac2 * rdm1->element(i2, i0);
        shalf->element(i1, i0, i2, i1) +=         rdm1->element(i2, i0);
      }
    }
  auto sview = group(group(*shalf, 2,4),0,2);
  shalf_hh_->copy_block(dim*jst, dim*ist, dim, dim, sview);

  shared_ptr<RDM<3,DataType>> r3 = rdm3->clone();
  sort_indices<1,3,0,2,4,0,1,1,1>(rdm3->data(), r3->data(), nact, nact, nact, nact, nact*nact);
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i4 = 0; i4 != nact; ++i4) {
        r3->element(i4, i2, i4, i3, i2, i3) +=  fac4 * fac;
        r3->element(i4, i2, i3, i4, i2, i3) += -fac2 * fac;
        r3->element(i4, i2, i2, i3, i4, i3) += -fac2 * fac;
        r3->element(i4, i2, i3, i2, i4, i3) +=  fac4 * fac;
        for (int i1 = 0; i1 != nact; ++i1) {
          r3->element(i1, i4, i4, i3, i2, i3) +=         rdm1->element(i2, i1);
          r3->element(i1, i4, i3, i4, i2, i3) += -fac2 * rdm1->element(i2, i1);
          r3->element(i1, i4, i3, i2, i4, i3) +=         rdm1->element(i2, i1);
          r3->element(i1, i4, i2, i3, i4, i3) += -fac2 * rdm1->element(i2, i1);
          r3->element(i1, i4, i1, i3, i2, i3) += -fac2 * rdm1->element(i2, i4);
          r3->element(i1, i2, i1, i4, i2, i3) += -fac2 * rdm1->element(i4, i3);
          r3->element(i4, i1, i3, i4, i2, i3) +=         rdm1->element(i2, i1);
          r3->element(i4, i1, i2, i4, i1, i3) +=         rdm1->element(i2, i3);
          r3->element(i2, i1, i3, i4, i2, i3) += -fac2 * rdm1->element(i4, i1);
          r3->element(i2, i1, i4, i3, i2, i3) +=         rdm1->element(i4, i1);
          r3->element(i2, i1, i1, i4, i2, i3) +=         rdm1->element(i4, i3);
          r3->element(i2, i1, i4, i1, i2, i3) += -fac2 * rdm1->element(i4, i3);
          r3->element(i4, i1, i4, i1, i2, i3) +=  fac4 * rdm1->element(i2, i3);
          r3->element(i4, i1, i1, i4, i2, i3) += -fac2 * rdm1->element(i2, i3);
          for (int i0 = 0; i0 != nact; ++i0) {
            r3->element(i0, i1, i3, i4, i2, i3) +=         rdm2->element(i2, i0, i4, i1);
            r3->element(i0, i1, i4, i3, i2, i3) +=         rdm2->element(i4, i0, i2, i1);
            r3->element(i0, i1, i1, i4, i2, i3) +=         rdm2->element(i4, i0, i2, i3);
            r3->element(i0, i1, i4, i1, i2, i3) += -fac2 * rdm2->element(i4, i0, i2, i3);
            r3->element(i0, i1, i2, i4, i1, i3) +=         rdm2->element(i2, i0, i4, i3);
            r3->element(i1, i0, i1, i4, i2, i3) += -fac2 * rdm2->element(i4, i0, i2, i3);
            r3->element(i1, i0, i4, i1, i2, i3) +=         rdm2->element(i4, i0, i2, i3);
            r3->element(i1, i0, i2, i4, i1, i3) +=         rdm2->element(i4, i0, i2, i3);
          }
        }
      }

  auto work2 = make_shared<MatType>(dim, dim);
  auto r3v = group(group(*r3, 4,6),0,4);
  auto workv = group(*work2, 0,2);
  contract(1.0, r3v, {0,1}, group(*fock_,0,2), {1}, 0.0, workv, {0});
  work_hh_->copy_block(dim*jst, dim*ist, dim, dim, work2);
}



template<typename DataType>
void Denom<DataType>::init_xh_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                             shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  auto shalf = make_shared<MatType>(dim*2, dim*2);
  RDM<2,DataType> ovl1 = *rdm2;
  RDM<2,DataType> ovl4 = *rdm2;
  ovl1.scale(-1.0);

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        ovl1.element(i, i, k, j) += fac2 * rdm1->element(k, j);
        ovl4.element(k, i, i, j) +=        rdm1->element(k, j);
      }

  MatType work(dim, dim);
  sort_indices<2,1,3,0,0,1,1,1>(ovl1.data(), work.data(), nact, nact, nact, nact);

  if (is_same<DataType,double>::value) {
    shalf->add_block(1.0, dim, dim, dim, dim, work);

    sort_indices<0,1,3,2,0,1,1,1>(ovl4.data(), work.data(), nact, nact, nact, nact);
    shalf->add_block(fac2, 0, 0, dim, dim, work);

    shalf->add_block(-1.0, dim, 0, dim, dim, work);
    shalf->add_block(-1.0, 0, dim, dim, dim, work);

    shalf_xh_->copy_block(dim*jst*2, dim*ist*2, dim*2, dim*2, shalf);
  } else {
    shalf_xh_->copy_block(dim*jst, dim*ist, dim, dim, work);
  }

  shared_ptr<RDM<3,DataType>> d0 = rdm3->copy();
  shared_ptr<RDM<3,DataType>> d3 = is_same<DataType,double>::value ? rdm3->copy() : nullptr;
  d0->scale(-1.0);
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4)
      for (int i3 = 0; i3 != nact; ++i3) {
        if (is_same<DataType,double>::value) {
          blas::ax_plus_y_n(1.0, rdm2->element_ptr(0,0, i4, i5), nact*nact, d3->element_ptr(0,0, i3, i5, i4, i3));
          blas::ax_plus_y_n(1.0, rdm1->element_ptr(0, i5), nact, d3->element_ptr(0, i4, i3, i5, i4, i3));
        }
        for (int i2 = 0; i2 != nact; ++i2) {
          if (is_same<DataType,double>::value) {
            blas::ax_plus_y_n(1.0, rdm2->element_ptr(0, i5, i2, i3), nact, d3->element_ptr(0, i4, i4, i5, i2, i3));
            blas::ax_plus_y_n(1.0, rdm2->element_ptr(0, i3, i4, i5), nact, d3->element_ptr(0, i2, i4, i5, i2, i3));
          }
          d0->element(i3, i2, i4, i5, i2, i3) +=  fac2 * rdm1->element(i4, i5);
          blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0, i3, i2, i5), nact, d0->element_ptr(0, i4, i2, i5, i4, i3));
          for (int i1 = 0; i1 != nact; ++i1) {
            d0->element(i3, i1, i4, i5, i2, i3) += -1.0  * rdm2->element(i2, i1, i4, i5);
            d0->element(i1, i1, i4, i5, i2, i3) +=  fac2 * rdm2->element(i2, i3, i4, i5);
          }
        }
      }

  MatType work2(dim, dim);

  auto d0v = group(group(*d0, 4,6),0,4);
  auto work2v = group(work2, 0,2);
  contract(1.0, d0v, {0,1}, group(*fock_,0,2), {1}, 0.0, work2v, {0});
  sort_indices<2,1,3,0,0,1,1,1>(work2.data(), work.data(), nact, nact, nact, nact);
  auto num = make_shared<MatType>(dim*2, dim*2);

  if (is_same<DataType,double>::value) {
    num->add_block(1.0, dim, dim, dim, dim, work);

    auto d3v = group(group(*d3, 4,6),0,4);
    contract(1.0, d3v, {0,1}, group(*fock_,0,2), {1}, 0.0, work2v, {0});
    sort_indices<0,1,3,2,0,1,1,1>(work2.data(), work.data(), nact, nact, nact, nact);
    num->add_block(fac2, 0, 0, dim, dim, work);

    num->add_block(-1.0, dim, 0, dim, dim, work);
    num->add_block(-1.0, 0, dim, dim, dim, work);

    work_xh_->copy_block(dim*jst*2, dim*ist*2, dim*2, dim*2, num);
  } else {
    work_xh_->copy_block(dim*jst, dim*ist, dim, dim, work);
  }
}


template<typename DataType>
void Denom<DataType>::init_xxh_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                              shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact*nact;
  shared_ptr<RDM<3,DataType>> ovl = rdm3->copy();
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i4 = 0; i4 != nact; ++i4)
        for (int i3 = 0; i3 != nact; ++i3)
          blas::ax_plus_y_n(1.0, rdm2->element_ptr(0, i4, i0, i5), nact, ovl->element_ptr(0, i3, i3, i4, i0, i5));
  auto shalf = make_shared<MatType>(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(ovl->data(), shalf->data(), nact, nact, nact, nact, nact, nact);
  shalf_xxh_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  // TODO a littile of duplication. Maybe should merge with xhh below
  shared_ptr<RDM<3,DataType>> fr4 = rdm4f->copy();
  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3) {
      const DataType f = fock_->element(i3, i4);
      for (int i7 = 0; i7 != nact; ++i7)
        for (int i0 = 0; i0 != nact; ++i0)
          for (int i6 = 0; i6 != nact; ++i6)
            blas::ax_plus_y_n(f, rdm2->element_ptr(0, i6, i0, i7), nact, fr4->element_ptr(0, i3, i4, i6, i0, i7));
    }
  // terms with frdm3(....) = rdm3(....ij) * f(i,j)
  {
    shared_ptr<RDM<2,DataType>> frdm3 = rdm2->clone();
    auto rdm3v = group(group(*rdm3,4,6),0,4);
    auto frdm3v = group(*frdm3,0,4);
    btas::contract(1.0, rdm3v, {0,1}, group(*fock_,0,2), {1}, 0.0, frdm3v, {0});
    for (int i7 = 0; i7 != nact; ++i7)
      for (int i0 = 0; i0 != nact; ++i0)
        for (int i6 = 0; i6 != nact; ++i6)
          for (int i5 = 0; i5 != nact; ++i5)
            blas::ax_plus_y_n(1.0, frdm3->element_ptr(0, i6, i0, i7), nact, fr4->element_ptr(0, i5, i5, i6, i0, i7));
  }
  // terms with grdm3(.....i) = rdm3(.....j) * f(i,j)
  {
    shared_ptr<RDM<3,DataType>> grdm3 = rdm3->clone();
    auto grdm3v = group(*grdm3,0,5);
    btas::contract(1.0, group(*rdm3,0,5), {0,1}, *fock_, {2,1}, 0.0, grdm3v, {0,2});
    sort_indices<1,0,1,1,1,1>(grdm3->data(), fr4->data(), nact*nact*nact*nact, nact*nact);
  }
  // terms with hrdm3(j.....) = rdm3(i.....) * f(i,j)
  {
    shared_ptr<RDM<3,DataType>> hrdm3 = rdm3->clone();
    auto hrdm3v = group(*hrdm3,1,6);
    btas::contract(1.0, group(*rdm3,1,6), {0,1}, *fock_, {0,2}, 0.0, hrdm3v, {2,1});
    sort_indices<1,0,2,1,1,1,1>(hrdm3->data(), fr4->data(), nact*nact, nact*nact, nact*nact);
  }
  auto fss = make_shared<MatType>(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(fr4->data(), fss->data(), nact, nact, nact, nact, nact, nact);
  work_xxh_->copy_block(dim*jst, dim*ist, dim, dim, fss);
}


template<typename DataType>
void Denom<DataType>::init_xhh_(const int jst, const int ist, shared_ptr<const RDM<1,DataType>> rdm1, shared_ptr<const RDM<2,DataType>> rdm2,
                                                              shared_ptr<const RDM<3,DataType>> rdm3, shared_ptr<const RDM<3,DataType>> rdm4f) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact*nact;
  shared_ptr<RDM<3,DataType>> ovl = rdm3->copy();
  ovl->scale(-1.0);
  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4) {
      blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0,0,0, i5), nact*nact*nact, ovl->element_ptr(0,0,0, i4, i4, i5));
      for (int i2 = 0; i2 != nact; ++i2) {
        blas::ax_plus_y_n(fac2, rdm2->element_ptr(0,0, i4, i5), nact*nact, ovl->element_ptr(0,0, i2, i2, i4, i5));
        blas::ax_plus_y_n(fac2, rdm1->element_ptr(0, i5), nact, ovl->element_ptr(0, i4, i2, i2, i4, i5));
        blas::ax_plus_y_n(-1.0, rdm1->element_ptr(0, i5), nact, ovl->element_ptr(0, i2, i2, i4, i4, i5));
        for (int i3 = 0; i3 != nact; ++i3) {
          blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0, i5, i3, i2), nact, ovl->element_ptr(0, i4, i3, i2, i4, i5));
          blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0, i2, i4, i5), nact, ovl->element_ptr(0, i3, i3, i2, i4, i5));
        }
      }
    }
  auto shalf = make_shared<MatType>(dim, dim);
  sort_indices<0,2,4,3,1,0,1,1,1>(ovl->data(), shalf->data(), nact*nact, nact, nact, nact, nact);
  shalf_xhh_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<3,DataType>> fr4 = rdm4f->copy();
  fr4->scale(-1.0);

  // terms with 1 or 2RDM intermediates
  {
    shared_ptr<RDM<2,DataType>> frdm3 = rdm2->clone();
    auto rdm3v = group(group(*rdm3,4,6),0,4);
    auto frdm3v = group(*frdm3,0,4);
    btas::contract(1.0, rdm3v, {0,1}, group(*fock_,0,2), {1}, 0.0, frdm3v, {0});

    shared_ptr<RDM<1,DataType>> frdm2 = rdm1->clone();
    auto rdm2v = group(group(*rdm2,2,4),0,2);
    auto frdm2v = group(*frdm2,0,2);
    btas::contract(1.0, rdm2v, {0,1}, group(*fock_,0,2), {1}, 0.0, frdm2v, {0});

    shared_ptr<RDM<2,DataType>> grdm2 = rdm2->clone();
    shared_ptr<RDM<2,DataType>> grdm2t = rdm2->clone();
    auto grdm2v = group(*grdm2,0,3);
    btas::contract(1.0, group(*rdm2,0,3), {0,1}, *fock_, {2,1}, 0.0, grdm2v, {0,2});
    sort_indices<1,0,0,1,1,1>(grdm2->data(), grdm2t->data(), nact*nact, nact*nact);

    shared_ptr<RDM<2,DataType>> hrdm2 = rdm2->clone();
    shared_ptr<RDM<2,DataType>> hrdm2t = rdm2->clone();
    auto hrdm2v = group(*hrdm2,1,4);
    btas::contract(1.0, group(*rdm2,1,4), {0,1}, *fock_, {0,2}, 0.0, hrdm2v, {2,1});
    sort_indices<1,0,0,1,1,1>(hrdm2->data(), hrdm2t->data(), nact*nact, nact*nact);

    for (int i7 = 0; i7 != nact; ++i7) {
      for (int i6 = 0; i6 != nact; ++i6) {
        blas::ax_plus_y_n(-1.0, grdm2t->element_ptr(0,0,0, i7), nact*nact*nact, fr4->element_ptr(0,0,0, i6, i6, i7));
        blas::ax_plus_y_n(-1.0, hrdm2t->element_ptr(0,0,0, i7), nact*nact*nact, fr4->element_ptr(0,0,0, i6, i6, i7));
        for (int i2 = 0; i2 != nact; ++i2) {
          blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0,0, i6, i7), nact*nact, fr4->element_ptr(0,0, i6, i2, i2, i7));
          blas::ax_plus_y_n(fac2, frdm3->element_ptr(0,0, i6, i7), nact*nact, fr4->element_ptr(0,0, i2, i2, i6, i7));
          blas::ax_plus_y_n(fac2, frdm2->element_ptr(0, i7), nact, fr4->element_ptr(0, i6, i2, i2, i6, i7));
          blas::ax_plus_y_n(-1.0, frdm2->element_ptr(0, i7), nact, fr4->element_ptr(0, i2, i2, i6, i6, i7));
          blas::ax_plus_y_n(fac2, grdm2t->element_ptr(0,0, i6, i7), nact*nact, fr4->element_ptr(0,0, i2, i2, i6, i7));
          blas::ax_plus_y_n(fac2, hrdm2t->element_ptr(0,0, i6, i7), nact*nact, fr4->element_ptr(0,0, i2, i2, i6, i7));
          for (int i5 = 0; i5 != nact; ++i5) {
            blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0, i7, i5, i2), nact, fr4->element_ptr(0, i6, i5, i2, i6, i7));
            blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0, i2, i6, i7), nact, fr4->element_ptr(0, i5, i5, i2, i6, i7));
            blas::ax_plus_y_n(-1.0,  grdm2->element_ptr(0, i7, i2, i5), nact, fr4->element_ptr(0, i6, i2, i5, i6, i7));
            blas::ax_plus_y_n(-1.0, grdm2t->element_ptr(0, i5, i6, i7), nact, fr4->element_ptr(0, i2, i2, i5, i6, i7));
            blas::ax_plus_y_n(-1.0, hrdm2t->element_ptr(0, i7, i5, i2), nact, fr4->element_ptr(0, i6, i5, i2, i6, i7));
            blas::ax_plus_y_n(-1.0, hrdm2t->element_ptr(0, i2, i5, i7), nact, fr4->element_ptr(0, i6, i6, i2, i5, i7));
          }
        }
      }
    }
  }
  // terms with grdm3(.....i) = rdm3(.....j) * f(i,j)
  {
    shared_ptr<RDM<3,DataType>> grdm3 = rdm3->clone();
    auto grdm3v = group(*grdm3,0,5);
    btas::contract(1.0, group(*rdm3,0,5), {0,1}, *fock_, {2,1}, 0.0, grdm3v, {0,2});
    sort_indices<1,0,1,1,-1,1>(grdm3->data(), fr4->data(), nact*nact*nact*nact, nact*nact);
    sort_indices<0,2,1,1,1,-1,1>(grdm3->data(), fr4->data(), nact*nact, nact*nact, nact*nact);
  }
  // terms with hrdm3(j.....) = rdm3(i.....) * f(i,j)
  {
    shared_ptr<RDM<3,DataType>> hrdm3 = rdm3->clone();
    auto hrdm3v = group(*hrdm3,1,6);
    btas::contract(1.0, group(*rdm3,1,6), {0,1}, *fock_, {0,2}, 0.0, hrdm3v, {2,1});
    sort_indices<1,0,1,1,-1,1>(hrdm3->data(), fr4->data(), nact*nact, nact*nact*nact*nact);
    sort_indices<1,0,2,1,1,-1,1>(hrdm3->data(), fr4->data(), nact*nact, nact*nact, nact*nact);
  }

  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3) {
      const DataType f = fock_->element(i3, i4);
      const DataType f2 = f * fac2;
      for (int i7 = 0; i7 != nact; ++i7) {
        blas::ax_plus_y_n(-f, rdm2->element_ptr(0,0,0, i7), nact*nact*nact, fr4->element_ptr(0,0,0, i3, i4, i7));
        for (int i6 = 0; i6 != nact; ++i6) {
          blas::ax_plus_y_n(f2, rdm2->element_ptr(0,0, i6, i7), nact*nact, fr4->element_ptr(0,0, i4, i3, i6, i7));
          blas::ax_plus_y_n(f2, rdm1->element_ptr(0, i7), nact, fr4->element_ptr(0, i6, i4, i3, i6, i7));
          blas::ax_plus_y_n(-f, rdm1->element_ptr(0, i7), nact, fr4->element_ptr(0, i6, i6, i3, i4, i7));
          blas::ax_plus_y_n(-f, rdm1->element_ptr(0, i7), nact, fr4->element_ptr(0, i3, i4, i6, i6, i7));
          blas::ax_plus_y_n(f2, rdm1->element_ptr(0, i7), nact, fr4->element_ptr(0, i3, i6, i6, i4, i7));
          for (int i2 = 0; i2 != nact; ++i2) {
            blas::ax_plus_y_n(-f, rdm2->element_ptr(0, i7, i6, i2), nact, fr4->element_ptr(0, i3, i6, i2, i4, i7));
            blas::ax_plus_y_n(-f, rdm2->element_ptr(0, i2, i6, i7), nact, fr4->element_ptr(0, i3, i4, i2, i6, i7));
          }
        }
      }
    }
  auto fss = make_shared<MatType>(dim, dim);
  sort_indices<0,2,4,3,1,0,1,1,1>(fr4->data(), fss->data(), nact*nact, nact, nact, nact, nact);
  work_xhh_->copy_block(dim*jst, dim*ist, dim, dim, fss);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class Denom<double>;
template class Denom<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
