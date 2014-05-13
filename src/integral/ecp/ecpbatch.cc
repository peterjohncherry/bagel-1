//
// BAGEL - Parallel electron correlation program.
// Filename: ecpbatch.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/integral/ecp/ecpbatch.h>
#include <src/integral/carsphlist.h>
#include <src/integral/sortlist.h>


using namespace bagel;
using namespace std;

const static CarSphList carsphlist;

ECPBatch::ECPBatch(const array<shared_ptr<const Shell>,2>& info, const shared_ptr<const Molecule> mol,
                         shared_ptr<StackMem> stack)
 : basisinfo_(info), mol_(mol) {

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

  integral_thresh_ = PRIM_SCREEN_THRESH;
  max_iter_ = 20;

  spherical_ = basisinfo_[0]->spherical();
  assert(spherical_ == basisinfo_[1]->spherical());

  common_init();

}


ECPBatch::~ECPBatch() {
  stack_->release(size_alloc_, stack_save_);
  if (allocated_here_) resources__->release(stack_);
}


void ECPBatch::compute() {

  const double zero = 0.0;
  const SortList sort(spherical_);

  fill_n(data_, size_alloc_, zero);
  double* const intermediate_c = stack_->get(cont0_ * cont1_ * asize_);
  double* current_data = intermediate_c;


  for (int contA = 0; contA != basisinfo_[0]->contractions().size(); ++contA) {
    for (int contC = 0; contC != basisinfo_[1]->contractions().size(); ++contC) {
      for (int izA = 0; izA <= ang0_; ++izA) {
        for (int iyA = 0; iyA <= ang0_ - izA; ++iyA) {
          const int ixA = ang0_ - izA - iyA;
          const array<int, 3> lA = {ixA, iyA, izA};
          for (int izC = 0; izC <= ang1_; ++izC) {
            for (int iyC = 0; iyC <= ang1_ - izC; ++iyC) {
              const int ixC = ang1_ - izC - iyC;
              const array<int, 3> lC = {ixC, iyC, izC};
              double tmp = 0.0;
              for (auto& aiter : mol_->atoms()) {
                shared_ptr<const ECP> aiter_ecp = aiter->ecp_parameters();

                RadialInt<AngularBatch, const shared_ptr<const ECP>, const array<shared_ptr<const Shell>,2>&,
                                        const double, const double, const array<int, 3>, const array<int, 3>>
                              radint(aiter_ecp, basisinfo_, contA, contC, lA, lC, false, max_iter_, integral_thresh_);
                tmp += radint.integral();
              }
              *current_data++ = tmp;
            }
          }
        }
      }
    }
  }

  double* const intermediate_fi = stack_->get(cont0_ * cont1_ * asize_intermediate_);
  const unsigned int array_size = cont0_ * cont1_ * asize_intermediate_;
  copy_n(intermediate_c, array_size, intermediate_fi);

  if (spherical_) {
    double* const intermediate_i = stack_->get(cont0_ * cont1_ * asize_final_);
    const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_fi, intermediate_i);

    const static SortList sort(true);
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data_, intermediate_i, cont1_, cont0_, 1, swap01_);
    stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
  } else {
    const static SortList sort(false);
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data_, intermediate_fi, cont1_, cont0_, 1, swap01_);
  }

  stack_->release(cont0_*cont1_*asize_intermediate_, intermediate_fi);
  stack_->release(cont0_*cont1_*asize_, intermediate_c);

}

void ECPBatch::common_init() {

  assert(basisinfo_.size() == 2);

  ang0_ = basisinfo_[0]->angular_number();
  ang1_ = basisinfo_[1]->angular_number();

  if (ang0_ < ang1_) {
    std::swap(basisinfo_[0], basisinfo_[1]);
    std::swap(ang0_, ang1_);
    swap01_ = true;
  } else {
    swap01_ = false;
  }

  cont0_ = basisinfo_[0]->num_contracted();
  cont1_ = basisinfo_[1]->num_contracted();

  amax_ = ang0_ + ang1_;
  amax1_ = amax_ + 1;
  amin_ = ang0_;

  asize_ = 0;
  for (int i = amin_; i != amax1_; ++i) asize_ += (i+1)*(i+2) / 2;
  asize_intermediate_ = (ang0_+1) * (ang0_+2) * (ang1_+1) * (ang1_+2) / 4;
  asize_final_ = (basisinfo_[0]->spherical() ? (2*ang0_+1) : (ang0_+1)*(ang0_+2)/2)
               * (basisinfo_[1]->spherical() ? (2*ang1_+1) : (ang1_+1)*(ang1_+2)/2);

  size_alloc_ = cont0_ * cont1_ * max(asize_intermediate_, asize_);

  stack_save_ = stack_->get(size_alloc_);
  data_ = stack_save_;

}
