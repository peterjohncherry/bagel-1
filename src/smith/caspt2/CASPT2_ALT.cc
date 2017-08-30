//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
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
#include <src/smith/caspt2/CASPT2_ALT.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in ) { 
////////////////////////////////////////////////////////////////////
  orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
}
////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(std::shared_ptr<const SMITH_Info<double>> ref){
////////////////////////////////////////////////////////////////////
  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  //det_ = ref->ciwfn()->det();
  //cc_ = ref->ciwfn()->civectors();
}

/////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::test() { 
/////////////////////////////////
  auto weqn = make_shared<Equation<Tensor>>();
  weqn->equation_build();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma2  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::compute_gamma12(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//  if (det_->compress()) { // uncompressing determinants
//    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
//    cc_->set_det(detex);
//  }
//  shared_ptr<Civec> ccbra = cc_->data(MM);
//  shared_ptr<Civec> ccket = cc_->data(NN);
 
 // shared_ptr<RDM<1>> gamma1;
 // shared_ptr<RDM<2>> gamma2;
 // tie(gamma1, gamma2) = compute_gamma12_from_civec(ccbra, ccket);

 // gamma1_->emplace(MM, NN, gamma1);
 // gamma2_->emplace(MM, NN, gamma2);
 //
 // cc_->set_det(det_); // 
 //
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
void
CASPT2_ALT::CASPT2_ALT::compute_gamma12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // since we consider here number conserving operators...
 // auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
 // sigma_2a1(cbra, dbra);
 // sigma_2a2(cbra, dbra);
 //
 // shared_ptr<Dvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
 // if (cbra != cket) {
 //   dket = make_shared<Dvec>(cket->det(), norb_*norb_);
 //   sigma_2a1(cket, dket);
 ////   sigma_2a2(cket, dket);
 //// } else {
 ////   dket = dbra;
 //// }
 ////
 // return compute_gamma12_last_step(dbra, dket, cbra);
}//

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
void
CASPT2_ALT::CASPT2_ALT::compute_gamma12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // const int nri = cibra->asize()*cibra->lenb();
 // const int ij  = norb_*norb_;
 //
 // // 1gamma c^dagger <I|\hat{E}|0>
 // // 2gamma \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
 // auto gamma1 = make_shared<RDM<1>>(norb_);
 // auto gamma2 = make_shared<RDM<2>>(norb_);
 // {
 //   auto cibra_data = make_shared<VectorB>(nri);
 //   copy_n(cibra->data(), nri, cibra_data->data());

 //   auto dket_data = make_shared<Matrix>(nri, ij);
 //   for (int i = 0; i != ij; ++i)
 //     copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
 //   auto gamma1t = btas::group(*gamma1,0,2);
 //   btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, gamma1t, {1});

 //   auto dbra_data = dket_data;
 //   if (dbra != dket) {
 //     dbra_data = make_shared<Matrix>(nri, ij);
 //     for (int i = 0; i != ij; ++i)
 //       copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
 //   }
 //   auto gamma2t = group(group(*gamma2, 2,4), 0,2);
 //   btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, gamma2t, {0,2});
 // }
 //
  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
 // unique_ptr<double[]> buf(new double[norb_*norb_]);
 // for (int i = 0; i != norb_; ++i) {
 //   for (int k = 0; k != norb_; ++k) {
 ////     copy_n(&gamma2->element(0,0,k,i), norb_*norb_, buf.get());
 ////     blas::transpose(buf.get(), norb_, norb_, gamma2->element_ptr(0,0,k,i));
 //   }
 // }
 //
 // return tie(gamma1, gamma2);
 return;
}
////////////////////////////////////////////////////////////////////
WICKTOOLS::WICKTOOLS::WICKTOOLS(std::shared_ptr<const SMITH_Info<double>> ref){
////////////////////////////////////////////////////////////////////
  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  // det_ = ref->ciwfn()->det();
  // cc_ = ref->ciwfn()->civectors();
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma2  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void WICKTOOLS::WICKTOOLS::compute_gamma12(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
void
WICKTOOLS::WICKTOOLS::compute_gamma12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
void
WICKTOOLS::WICKTOOLS::compute_gamma12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 return;
}
#endif
