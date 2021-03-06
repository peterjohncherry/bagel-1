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


#include <iostream>
#include <iomanip>
#include <src/smith/caspt2/CASPT2.h>
#include <src/util/math/linearRM.h>
#include <src/smith/caspt2/MSCASPT2.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/integrals/moint_computer.h>
#include <src/prop/proptool/debugging_utils.h>
#define __TEST_PROPTOOL
#ifdef  __TEST_PROPTOOL
#include <src/prop/proptool/proptool.h>
#endif

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = info_->nact() ? ref->ciwfn()->nstates() : 1;
  
  // MS-CASPT2: t2 and s as MultiTensor (t2all, sall)
  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (int j = 0; j != nstates_; ++j)
      if (!info_->sssr() || i == j) {
        (*tmp)[j] = init_amplitude();
        (*tmp2)[j] = init_residual();
      }
    t2all_.push_back(tmp);
    sall_.push_back(tmp2);
    rall_.push_back(tmp2->clone());
  }

  energy_.resize(nstates_);
  pt2energy_.resize(nstates_);
}


CASPT2::CASPT2::CASPT2(const CASPT2& cas) : SpinFreeMethod(cas) {
cout << "CASPT2::CASPT2::CASPT2(const CASPT2& cas) : SpinFreeMethod(cas) {" << endl;
  info_    = cas.info_;
  closed_  = cas.closed_;
  rvirt_   = cas.rvirt_;
  ractive_ = cas.ractive_;
  rclosed_ = cas.rclosed_;
  nstates_ = cas.nstates_;
  nstates_ = 1;// cas.nstates_;
  heff_    = cas.heff_;
  fockact_ = cas.fockact_;
  e0all_   = cas.e0all_;
  xmsmat_  = cas.xmsmat_;
  energy_  = cas.energy_;
  pt2energy_ = cas.pt2energy_;

  // sall is changed in gradient and nacme codes while the others are not
  t2all_ = cas.t2all_;
  rall_  = cas.rall_;
  for (int i = 0; i != nstates_; ++i) {
    sall_.push_back(cas.sall_[i]->copy());
  }
  h1_ = cas.h1_;
  f1_ = cas.f1_;
  v2_ = cas.v2_;
  H_2el_ =  cas.H_2el_; 

  rdm0all_ = cas.rdm0all_;
  rdm1all_ = cas.rdm1all_;
  rdm2all_ = cas.rdm2all_;
  rdm3all_ = cas.rdm3all_;
  rdm4all_ = cas.rdm4all_;

}

void CASPT2::CASPT2::do_rdm_deriv(double factor) {
  Timer timer(1);
  tie(den0cirdmt, den1cirdmt, den2cirdmt, den3cirdmt, den4cirdmt) = feed_denci();

  ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), 1);
  const size_t nact  = info_->nact();
  const size_t norb2 = nact * nact;
  const size_t ndet = ci_deriv_->data(0)->size();
  const size_t ijmax = info_->cimaxchunk();
  const size_t ijnum = ndet * norb2 * norb2;
  const size_t npass = ((mpi__->size() > ((ijnum - 1)/ijmax + 1)) && (mpi__->size() != 1) && ndet > 10000) ? mpi__->size() : (ijnum - 1) / ijmax + 1;
  const size_t nsize = (ndet - 1) / npass + 1;

  if (npass > 1)
    cout << "       - CI derivative contraction will be done with " << npass << " passes" << endl;

  // Fock-weighted 2RDM derivative evaluated first (needed for calculating Fock-weighted 3RDM derivative)
  rdm2fderiv_ = SpinFreeMethod<double>::feed_rdm_2fderiv(info_, fockact_, 0);

  if (npass > 1)
    timer.tick_print("Fock-weighted 2RDM derivative");

  // embarrasingly parallel mode. npass > 1 -> distribute among the nodes.
  // otherwise just do using all the nodes.
  const int nproc = npass > 1 ? 1 : mpi__->size();
  const int ncomm = mpi__->size() / nproc;
  const int icomm = mpi__->rank() / nproc;
  mpi__->split(nproc);

  for (int ipass = 0; ipass != npass; ++ipass) {
    if (ipass % ncomm == icomm && ncomm != icomm) {
      const size_t ioffset = ipass * nsize;
      const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;
      tie(rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3fderiv_)
        = SpinFreeMethod<double>::feed_rdm_deriv(info_, fockact_, 0, ioffset, isize, rdm2fderiv_);

      shared_ptr<VectorB> bdata = contract_rdm_deriv(info_->ciwfn(), ioffset, isize, fockact_);
      blas::ax_plus_y_n(factor, bdata->data(), ndet, ci_deriv_->data(0)->data());

      if (npass > 1) {
        stringstream ss; ss << "Multipassing (" << setw(2) << ipass + 1 << " / " << npass << ")";
        timer.tick_print(ss.str());
      }
    }
  }
  mpi__->merge();

  if (npass > 1)
    mpi__->allreduce(ci_deriv_->data(0)->data(), ndet);
}


void CASPT2::CASPT2::solve() {
cout << "CASPT2::CASPT2::solve" << endl;
  Timer timer;
  print_iteration();

  const int ncore = info_->ncore();
  const int nclosed = info_->nclosed()-info_->ncore();
  const int nact = info_->nact();
  const int nvirt = info_->nvirt() + info_->nfrozenvirt();
 
  // <proj_jst|H|0_K> set to sall in ms-caspt2
  cout << "getting  <proj_jst|H|0_K> " << endl;
  for (int istate = 0; istate != nstates_; ++istate) { //K states
    t2all_[istate]->fac(istate) = 0.0;
    sall_[istate]->fac(istate)  = 0.0;

    for (int jst=0; jst != nstates_; ++jst) { // <jst|
      if (info_->sssr() && jst != istate)
        continue;
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
      while(!sourceq->done())
        sourceq->next_compute();
    }
  }
  cout << "got  <proj_jst|H|0_K> " << endl;

  // solve linear equation for t amplitudes
  t2all_ = solve_linear(sall_, t2all_);
  timer.tick_print("CASPT2 energy evaluation");
  cout << endl;

  for (int istate = 0; istate != nstates_; ++istate) {
    n = init_residual();
    double norm = 0.0;
    for (int jst = 0; jst != nstates_; ++jst) { // bra
      for (int ist = 0; ist != nstates_; ++ist) { // ket
        if (info_->sssr() && (jst != istate || ist != istate))
          continue;
        set_rdm(jst, ist);
        t2 = t2all_[istate]->at(ist);
        shared_ptr<Queue> normq = make_normq(true, jst == ist);
        while (!normq->done())
          normq->next_compute();
        norm += dot_product_transpose(n, t2all_[istate]->at(jst));
      }
    }

    pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate) - info_->shift()*norm;
    cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    if (info_->shift() != 0.0)
      cout << "        w/o shift correction  " << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) << endl;
    cout << "        reference weight      " << fixed << setw(20) << setprecision(10) << 1.0/(1.0+norm) << endl;
    cout << endl;
  }

  // TODO Implement off-diagonal shift correction for nonrelativistic energy + nuclear gradients
  if (info_->shift() && info_->do_ms() && !info_->shift_diag())
    cout << "    Applying levelshift correction to diagonal elements of the Hamiltonian only.  (Off-diagonals have not been implemented for non-relativistic CASPT2.)" << endl << endl;

  // MS-CASPT2
  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<Matrix>(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      auto sist = make_shared<MultiTensor>(nstates_);
      for (int jst = 0; jst != nstates_; ++jst) {
        if (sall_[ist]->at(jst)) {
          sist->at(jst) = sall_[ist]->at(jst);
        } else {
          set_rdm(jst, ist);
          s = init_residual();
          shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
          while(!sourceq->done())
            sourceq->next_compute();
          sist->at(jst) = s;
        }
      }

      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          (*heff_)(ist, ist) = pt2energy_[ist];
        } else {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > ]
          (*heff_)(jst, ist) = dot_product_transpose(sist, t2all_[jst]) + (*eref_)(jst, ist);
        }
      }
    }
    heff_->symmetrize();

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASPT2 Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    heff_->diagonalize(eig);
    copy_n(eig.data(), nstates_, pt2energy_.data());

    // print out the eigen vector
    cout << endl;
    cout << "    * MS-CASPT2 rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    if (xmsmat_) {
      cout << endl;
      cout << "    * XMS-CASPT2 rotation matrix";
      for (int ist = 0; ist != nstates_; ++ist) {
        cout << endl << "      ";
        for (int jst = 0; jst != nstates_; ++jst)
          cout << setw(20) << setprecision(10) << msrot()->element(ist, jst);
      }
      cout << endl << endl;
    }

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate)
      cout << "    * MS-CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    cout << endl << endl;
  } else {
    heff_ = make_shared<Matrix>(1,1);
    heff_->element(0,0) = 1.0;
  }
  energy_ = pt2energy_;
   
  SMITH::IndexRange brng_act =  *ractive_;
  SMITH::IndexRange brng_core = *rclosed_;
  SMITH::IndexRange brng_virt = *rvirt_;
  
  SMITH::IndexRange brng_free = *rclosed_;
  brng_free.merge(*ractive_);
  brng_free.merge(*rvirt_);
 
#ifdef __TEST_PROPTOOL
  SMITH::IndexRange brng_occ = *rclosed_;
  brng_occ.merge(*ractive_);

  SMITH::IndexRange brng_exc = *ractive_;
  brng_exc.merge(*rvirt_);

   
  vector<SMITH::IndexRange> branges = { brng_core, brng_act, brng_virt, brng_free, brng_occ, brng_exc }; 
  vector<char> range_names = { 'c', 'a', 'v', 'f', 'o', 'e' }; 
  for ( int r0 = 0 ; r0!=6; r0++) { 
  for ( int r1 = 0 ; r1!=6; r1++) { 
  for ( int r2 = 0 ; r2!=6; r2++) { 
  for ( int r3 = 0 ; r3!=6; r3++) { 
    string rname = "";
    rname += range_names[r0];rname += range_names[r1]; rname += range_names[r2]; rname += range_names[r3];
    vector<SMITH::IndexRange> rblock = { branges[r0], branges[r1], branges[r2], branges[r3] }; 
    range_block_map_.emplace( rname, rblock );
  }}}}
  
  cout << endl <<  "======================  BEGIN PROPTOOL TEST =========================" << endl;
  shared_ptr<PropTool::PropTool> proptool = make_shared<PropTool::PropTool>( proptool_input_, info_->geom(), info_->ref() );
  cout << "S::cpt2::X1" << endl;
   
  shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> tensor_calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>();
  auto rcm = make_shared<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>>();
  auto closed_rng  = make_shared<SMITH::IndexRange>(*rclosed_);
  auto active_rng  = make_shared<SMITH::IndexRange>(*ractive_);
  auto virtual_rng = make_shared<SMITH::IndexRange>(*rvirt_);
  
  auto free_rng   = make_shared<SMITH::IndexRange>(*closed_rng);
  free_rng->merge(*active_rng);
  free_rng->merge(*virtual_rng);
  
  auto not_closed_rng  = make_shared<SMITH::IndexRange>(*active_rng); not_closed_rng->merge(*virtual_rng);
  auto not_active_rng  = make_shared<SMITH::IndexRange>(*closed_rng); not_active_rng->merge(*virtual_rng);
  auto not_virtual_rng = make_shared<SMITH::IndexRange>(*closed_rng); not_virtual_rng->merge(*active_rng);

  rcm->emplace("c", closed_rng);
  rcm->emplace("a", active_rng);
  rcm->emplace("v", virtual_rng);
  rcm->emplace("free", free_rng);
  rcm->emplace("notcor",not_closed_rng);
  rcm->emplace("notact",not_active_rng);
  rcm->emplace("notvir", not_virtual_rng);
  proptool->set_range_conversion_map( rcm );
  vector<SMITH::IndexRange> rng_block_nz = { *(rcm->at("c")), *(rcm->at("v")), *(rcm->at("c")), *(rcm->at("v")) };

  tensor_calc->set_tensor_elems(t2all_[0]->at(0), 1.0);
  shared_ptr<SMITH::Tensor_<double>> tamps_full = build_full_T( t2all_[0]->at(0) ); 
  tensor_calc->zero_all_but_block( tamps_full ,  rng_block_nz );

  tensor_calc->set_tensor_elems(v2_, 1.0);
  tensor_calc->zero_all_but_block( v2_ ,  rng_block_nz );

  proptool->tamps_smith_ = tamps_full->copy();
  proptool->set_maxtile( info_->maxtile() ); 

  proptool->construct_task_lists();

  cout << "t2all_[0]->at(0)->norm() = " <<  t2all_[0]->at(0)->norm() << endl;
  cout << "tamps_full->norm() = " <<  tamps_full->norm() << endl;
  cout << "SMITH v2_->norm() = " <<  v2_->norm() << endl << endl << endl;
  
  cout << endl << endl <<endl;
  proptool->execute_compute_lists();

  cout << " xxxxxxx FROM moint_computer xxxxxxxxxx" << endl;
  cout << "            v2_->norm() = " << v2_->norm() << endl <<endl; 

  {// TEST source
 
    set_rdm(0, 0);
    double source_norm = 0.0;
    s = init_residual();
    s->zero();
    h1_->zero();
    shared_ptr<Queue> source_task_list = make_sourceq(true, true);
    while(!source_task_list->done())
      source_task_list->next_compute();

    cout << "v2 ranges = [ " ; cout.flush(); for (auto elem : v2_->indexrange() ) {cout << elem.size() << " " ; cout.flush(); } cout << " ] " << endl; 
    cout << "s  ranges = [ " ; cout.flush(); for (auto elem : s->indexrange() ) {cout << elem.size() << " " ; cout.flush(); } cout << " ] " << endl; 
    cout << "----------------------------------TEST SMITH-------------------------------" << endl;
    cout <<" dot_product_transpose(s, t2_one     ) = " <<  dot_product_transpose(s, t2all_[0]->at(0))<< endl; // + (*eref_)(0, 0);
    cout <<" dot_product_transpose(s, tamps_full ) = " <<  dot_product_transpose(s, tamps_full )<< endl; // + (*eref_)(0, 0);
    cout <<" dot_product(s, t2_one     )           = " <<  s->dot_product( t2all_[0]->at(0) ) << endl; // + (*eref_)(0, 0);
    cout <<" dot_product(s, tamps_full )           = " <<  s->dot_product( tamps_full ) << endl; // + (*eref_)(0, 0);
    cout << "post dot source_norm = "<<  s->norm() << endl;
    cout << "---------------------------------------------------------------------------" << endl;
  } //END TEST
  throw logic_error( "die here for testing purposes!" ); 
  #endif
}


// function to solve linear equation
vector<shared_ptr<MultiTensor_<double>>> CASPT2::CASPT2::solve_linear(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  Timer mtimer;
  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;
  for (int i = 0; i != nstates_; ++i) {  // K states
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i] - info_->shift();
    energy_[i] = 0.0;
    // set guess vector
    t[i]->zero();
    if (s[i]->rms() < 1.0e-15) {
      print_iteration(0, 0.0, 0.0, mtimer.tick());
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      update_amplitude(t[i], s[i]);
    }

    auto solver = make_shared<LinearRM<MultiTensor>>(info_->davidson_subspace(), s[i]);
    for (int iter = 0; iter != info_->maxiter(); ++iter) {
      rall_[i]->zero();

      const double norm = t[i]->norm();
      t[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
        for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
          if (info_->sssr() && (jst != i || ist != i))
            continue;
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t[i]->at(ist);
          r = rall_[i]->at(jst);

          // compute residuals named r for each K
          e0_ = e0all_[i] - info_->shift();
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
        }
      }
      // solve using subspace updates
      rall_[i] = solver->compute_residual(t[i], rall_[i]);
      t[i] = solver->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error = rall_[i]->norm() / pow(rall_[i]->size(), 0.25);
      print_iteration(iter, energy_[i], error, mtimer.tick());
      conv = error < info_->thresh();

      // compute delta t2 and update amplitude
      if (!conv) {
        t[i]->zero();
        update_amplitude(t[i], rall_[i]);
      }
      if (conv) break;
    }
    if (i+1 != nstates_) cout << endl;
    converged &= conv;
  }
  print_iteration(!converged);
  return t;
}


void CASPT2::CASPT2::solve_dm(const int istate, const int jstate) {
  {
    MSCASPT2::MSCASPT2 ms(*this);
    ms.solve_dm(istate, jstate);
    vden1_ = ms.vden1();
  }
}


void CASPT2::CASPT2::solve_gradient(const int targetJ, const int targetI, shared_ptr<const NacmType> nacmtype, const bool nocider) {
  Timer timer;
  // First solve lambda equation if this is MS-CASPT2
  assert (!((targetJ != targetI) && (nstates_ == 1)));

  if ((info_->do_ms() && nstates_ > 1) || info_->shift() != 0.0) {
    // Lambda equation solver
    for (int i = 0; i != nstates_; ++i)
      lall_.push_back(t2all_[i]->clone());
    print_iteration();

    // source stores the result of summation over M'
    if (targetJ == targetI) {
      // Gradient: is special case of targetJ = targetI.
      const int target = targetJ;

      auto source = make_shared<MultiTensor>(nstates_);
      for (auto& i : *source)
        i = init_residual();
      for (int ist = 0; ist != nstates_; ++ist) {//N states
        auto sist = make_shared<MultiTensor>(nstates_);
        for (int jst = 0; jst != nstates_; ++jst) {
          if (sall_[ist]->at(jst)) {
            sist->at(jst) = sall_[ist]->at(jst);
          } else {
            set_rdm(jst, ist);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
            while(!sourceq->done())
              sourceq->next_compute();
            sist->at(jst) = s;
          }
        }
        source->ax_plus_y((*heff_)(ist, target), sist);
      }

      for (int istate = 0; istate != nstates_; ++istate) { //L states
        sall_[istate]->zero();
        for (int jst = 0; jst != nstates_; ++jst)
          if (!info_->sssr() || istate == jst)
            sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, target), source->at(jst));
        if (info_->shift() != 0.0) {
          // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
          n = init_residual();
          for (int jst = 0; jst != nstates_; ++jst) { // bra
            for (int ist = 0; ist != nstates_; ++ist) { // ket
              if (info_->sssr() && (jst != istate || ist != istate))
                continue;
              set_rdm(jst, ist);
              t2 = t2all_[istate]->at(ist);
              shared_ptr<Queue> normq = make_normq(true, jst == ist);
              while (!normq->done())
                normq->next_compute();
              sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * pow((*heff_)(istate, target), 2.0), n);
            }
          }
        }
      }
    } else {
      // NACME case
      auto sourceJ = make_shared<MultiTensor>(nstates_);
      auto sourceI = make_shared<MultiTensor>(nstates_);
      for (auto& i : *sourceJ)
        i = init_residual();
      for (auto& i : *sourceI)
        i = init_residual();

      for (int ist = 0; ist != nstates_; ++ist) { // L states
        auto sist = make_shared<MultiTensor>(nstates_);
        for (int jst = 0; jst != nstates_; ++jst) {
          if (sall_[ist]->at(jst)) {
            sist->at(jst) = sall_[ist]->at(jst);
          } else {
            set_rdm(jst, ist);
            s = init_residual();
            shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
            while(!sourceq->done())
              sourceq->next_compute();
            sist->at(jst) = s;
          }
        }
        sourceJ->ax_plus_y((*heff_)(ist, targetI) * 0.5, sist);
        sourceI->ax_plus_y((*heff_)(ist, targetJ) * 0.5, sist);
      }

      for (int istate = 0; istate != nstates_; ++istate) { //K states
        sall_[istate]->zero();
        for (int jst = 0; jst != nstates_; ++jst)
          if (!info_->sssr() || istate == jst) {
            sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetI), sourceI->at(jst));
            sall_[istate]->at(jst)->ax_plus_y((*heff_)(istate, targetJ), sourceJ->at(jst));
          }
        if (info_->shift() != 0.0) {
          // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
          n = init_residual();
          for (int jst = 0; jst != nstates_; ++jst) { // bra
            for (int ist = 0; ist != nstates_; ++ist) { // ket
              if (info_->sssr() && (jst != istate || ist != istate))
                continue;
              set_rdm(jst, ist);
              t2 = t2all_[istate]->at(ist);
              shared_ptr<Queue> normq = make_normq(true, jst == ist);
              while (!normq->done())
                normq->next_compute();
              sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * (*heff_)(istate, targetJ) * (*heff_)(istate, targetI), n);
            }
          }
        }
      }
    }
    // solve linear equation and store lambda in lall
    lall_ = solve_linear(sall_, lall_);

    timer.tick_print("CASPT2 lambda equation");
  }

  if (lall_.empty()) {
    t2 = t2all_[0]->at(0);
    {
      den2 = h1_->clone();
      shared_ptr<Queue> dens2 = make_densityq();
      while (!dens2->done())
        dens2->next_compute();
      den2_ = den2->matrix();
    } {
      den1 = h1_->clone();
      shared_ptr<Queue> dens1 = make_density1q();
      while (!dens1->done())
        dens1->next_compute();
      den1_ = den1->matrix();
    } {
      Den1 = init_residual();
      shared_ptr<Queue> Dens1 = make_density2q();
      while (!Dens1->done())
        Dens1->next_compute();
      Den1_ = Den1;
    }
    timer.tick_print("Correlated density matrix evaluation");
 

    // then form deci0 - 4
    den0ci = rdm0_->clone();
    den1ci = rdm1_->clone();
    den2ci = rdm2_->clone();
    den3ci = rdm3_->clone();
    den4ci = rdm3_->clone();
    shared_ptr<Queue> dec = make_deciq(/*reset = */true);
    while (!dec->done())
      dec->next_compute();
    timer.tick_print("CI derivative evaluation");

    // when active is divided into the blocks, den4ci is evaluated (activeblock)**2 times
    double den4factor = 1.0 / static_cast<double>(active_.nblock() * active_.nblock());
    den4ci->scale(den4factor);

    // and contract them with rdm derivs
    do_rdm_deriv(1.0);

    timer.tick_print("CI derivative contraction");
    cout << endl;
  } else {
    // in case when CASPT2 is not variational...
    MSCASPT2::MSCASPT2 ms(*this);
    
    ms.solve_gradient(targetJ, targetI, nocider);
    den1_ = ms.rdm11();
    den2_ = ms.rdm12();
    Den1_ = ms.rdm21();
    if (!nocider)
      ci_deriv_ = ms.ci_deriv();
    dcheck_ = ms.dcheck();
    if (targetJ != targetI)
      vden1_ = ms.vden1();
    timer.tick();
  }

  correlated_norm_.resize(nstates_);
  if (nstates_ == 1 && info_->shift() == 0.0) {
    n = init_residual();
    shared_ptr<Queue> normq = make_normq();
    while (!normq->done())
      normq->next_compute();
    correlated_norm_[0] = dot_product_transpose(n, t2);
  } else {
    n = init_residual();
    for (int istate = 0; istate != nstates_; ++istate) {
      double tmp = 0.0;
      for (int jst = 0; jst != nstates_; ++jst) { // bra
        for (int ist = 0; ist != nstates_; ++ist) { // ket
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          shared_ptr<Queue> normq = make_normq(true, jst == ist);
          while (!normq->done())
            normq->next_compute();
          tmp += dot_product_transpose(n, lall_[istate]->at(jst));
        }
      }
      correlated_norm_[istate] = tmp;
    }
  }
  timer.tick_print("T1 norm evaluation");

  // some additional terms to be added
  const int ncore = info_->ncore();
  const int nclosed = info_->nclosed()-info_->ncore();
  const int nact = info_->nact();
  {
    // d_1^(2) -= <1|1><0|E_mn|0>     [Celani-Werner Eq. (A6)]
    auto dtmp = den2_->copy();
    for (int ist = 0; ist != nstates_; ++ist) {
      auto rdmtmp = rdm1all_->at(ist, ist)->matrix();
      for (int i = nclosed; i != nclosed+nact; ++i)
        for (int j = nclosed; j != nclosed+nact; ++j)
          dtmp->element(j, i) -=  correlated_norm_[ist] * (*rdmtmp)(j-nclosed, i-nclosed);
    }
    dtmp->symmetrize();
    den2_ = dtmp;
  }

  shared_ptr<const Reference> ref = info_->ref();
  const MatView acoeff = coeff_->slice(nclosed+ncore, nclosed+ncore+nact);

  // code to calculate h+g(d). When add is false, h is not added
  auto focksub = [&](shared_ptr<const Matrix> moden, const MatView coeff, const bool add) {
    shared_ptr<const Matrix> jop = ref->geom()->df()->compute_Jop(make_shared<Matrix>(coeff * *moden ^ coeff));
    auto out = make_shared<Matrix>(acoeff % (add ? (*ref->hcore() + *jop) : *jop) * acoeff);
    shared_ptr<const DFFullDist> full = ref->geom()->df()->compute_half_transform(acoeff)->apply_J()->compute_second_transform(coeff)->swap();
    shared_ptr<DFFullDist> full2 = full->copy();
    full2 = full2->transform_occ1(moden);
    *out += *full->form_2index(full2, -0.5);
    return out;
  };

  if (!nocider) {
    shared_ptr<const Matrix> fock = focksub(ref->rdm1_mat(), coeff_->slice(0, ref->nocc()), true); // f
    shared_ptr<const Matrix> gd2 = focksub(den2_, coeff_->slice(ncore, coeff_->mdim()), false); // g(d2)

    // correct cideriv for fock derivative [Celani-Werner Eq. (C1), some terms in first and second lines]
    // y_I += (g[d^(2)]_ij - Nf_ij) <I|E_ij|0>
    for (int ist = 0; ist != nstates_; ++ist) {
      const Matrix op(*gd2 * (1.0/nstates_) - *fock * correlated_norm_[ist]);
      shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
      for (int i = 0; i != nact; ++i)
        for (int j = 0; j != nact; ++j)
          ci_deriv_->data(ist)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
    }

    // y_I += <I|H|0> (for mixed states); taking advantage of the fact that unrotated CI vectors are eigenvectors
    if (targetJ == targetI) {
      // Gradient case. Special case of NACME with targetJ = targetI
      const Matrix ur(xmsmat_ ? *xmsmat_ * *heff_ : *heff_);
      const int target = targetJ;
      for (int ist = 0; ist != nstates_; ++ist)
        for (int jst = 0; jst != nstates_; ++jst)
          ci_deriv_->data(jst)->ax_plus_y(2.0*ur(ist,target)*(*heff_)(jst,target)*ref->energy(ist), info_orig_->ciwfn()->civectors()->data(ist));
    } else {
      // NACME case. targetJ and target I are separately used
      const Matrix ur(xmsmat_ ? *xmsmat_ * *heff_ : *heff_);
      for (int ist = 0; ist != nstates_; ++ist)
        for (int jst = 0; jst != nstates_; ++jst) {
          double urheff = (ur(ist,targetJ)*(*heff_)(jst,targetI) + ur(ist, targetI)*(*heff_)(jst,targetJ)) * ref->energy(ist);
          ci_deriv_->data(jst)->ax_plus_y(urheff, info_orig_->ciwfn()->civectors()->data(ist));
        }
    }

    // finally if this is XMS-CASPT2 gradient computation, we compute dcheck and contribution to y
    if (xmsmat_) {
      Matrix wmn(nstates_, nstates_);
      shared_ptr<Tensor> dc = rdm1_->clone();
      for (int i = 0; i != nstates_; ++i)
        for (int j = 0; j != i; ++j) {
          double cy = info_->ciwfn()->civectors()->data(j)->dot_product(ci_deriv_->data(i))
                    - info_->ciwfn()->civectors()->data(i)->dot_product(ci_deriv_->data(j));
          // If this is Full NACME, <U | dU/dX> contribution should be added
          if ((targetJ != targetI) && (nacmtype->is_full() || nacmtype->is_etf())) {
            cy += (pt2energy_[targetI] - pt2energy_[targetJ])
                * ((*heff_)(i,targetI) * (*heff_)(j,targetJ) - (*heff_)(j,targetI) * (*heff_)(i,targetJ));
          }
          wmn(j,i) = fabs(e0all_[j]-e0all_[i]) > 1.0e-12 ? -0.5 * cy / (e0all_[j]-e0all_[i]) : 0.0;
          wmn(i,j) = wmn(j,i);
          dc->ax_plus_y(wmn(j,i), rdm1all_->at(j, i));
          dc->ax_plus_y(wmn(i,j), rdm1all_->at(i, j));
        }
      dcheck_ = dc->matrix();

      // fill this into CI derivative. (Y contribution is done inside Z-CASSCF together with frozen core)
      shared_ptr<const Matrix> gdc = focksub(dcheck_, acoeff, false);
      for (int ist = 0; ist != nstates_; ++ist) {
        shared_ptr<const Dvec> deriv = ref->rdm1deriv(ist);
        for (int jst = 0; jst != nstates_; ++jst) {
          Matrix op(*fock * wmn(jst, ist));
          if (ist == jst)
            op += *gdc * (1.0/nstates_) * 0.5;
          for (int i = 0; i != nact; ++i)
            for (int j = 0; j != nact; ++j)
              ci_deriv_->data(jst)->ax_plus_y(2.0*op(j,i), deriv->data(j+i*nact));
        }
      }

      // also rotate cideriv back to the MS states
      btas::contract(1.0, *ci_deriv_->copy(), {0,1,2}, (*xmsmat_), {3,2}, 0.0, *ci_deriv_, {0,1,3});
    }
  }

  // restore original energy
  energy_ = pt2energy_;
  timer.tick_print("Postprocessing SMITH");
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TEST
shared_ptr<SMITH::Tensor_<double>> CASPT2::CASPT2::build_full_T( shared_ptr<SMITH::Tensor_<double>> smith_t) {
cout << " CASPT2::CASPT2::build_full_T( shared_ptr<SMITH::Tensor_<double>>) " << endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

   auto tensor_calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>(); 

   shared_ptr<SMITH::Tensor_<double>> t_amps_full = make_shared<SMITH::Tensor_<double>>(range_block_map_.at( "oeoe") );
   t_amps_full->allocate();
    
   vector<string> t_blocks = { "cvcv", "cvav", "avcv", "avav", "caca", "caaa", "aaca", "cacv", "caav", "aacv", "aaav", "cvca", "cvaa", "avca", "avaa" }; 

   vector<vector<int>> symm_orders = {{0,3,2,1},{2,1,0,3},{2,3,0,1}};  
   vector<double> symm_facs = { -1.0, -1,0, 1.0 };  
 
   auto reorder_vec = []( vector<SMITH::IndexRange> tbr, vector<int> ro ) { 
     vector<SMITH::IndexRange> nbr = { tbr[ro[0]] , tbr[ro[1]] , tbr[ro[2]], tbr[ro[3]] }; 
     return nbr; 
   }; 

   for ( auto block_name : t_blocks ) {
     cout<< "block_name = " << block_name << endl;
     vector<SMITH::IndexRange> t_block_ranges = range_block_map_.at(block_name);
     vector<SMITH::Index> t_block_index = { t_block_ranges[0].range(0), t_block_ranges[1].range(0), t_block_ranges[2].range(0), t_block_ranges[3].range(0) };

     if (smith_t->exists(t_block_index)) {
       shared_ptr<SMITH::Tensor_<double>> t_part = Tensor_Arithmetic_Utils::get_sub_tensor( smith_t, t_block_ranges );
       tensor_calc->put_sub_tensor( t_part, t_amps_full );

     } else {
       cout << " block [ "; cout.flush(); for( auto idx : t_block_ranges) { cout << idx.size() << " "; cout.flush(); } cout << "] not in smith_t" << endl;

       for ( int ii = 0; ii != 3; ++ii ) {

         vector<SMITH::IndexRange> reord_t_block_ranges = reorder_vec( t_block_ranges, symm_orders[ii] ); 
         cout << "trying block [ "; cout.flush(); for( auto idx : t_block_ranges) { cout << idx.size() << " "; cout.flush(); } cout << "] not in smith_t" << endl;
         vector<SMITH::Index> reord_t_block_index = { reord_t_block_ranges[0].range(0), reord_t_block_ranges[1].range(0), reord_t_block_ranges[2].range(0), reord_t_block_ranges[3].range(0) };

         if ( smith_t->exists(reord_t_block_index) ) { 
           cout << "Found!!" << endl;
           shared_ptr<SMITH::Tensor_<double>> tmp = Tensor_Arithmetic_Utils::get_sub_tensor( smith_t, reord_t_block_ranges );
           auto t_part = tensor_calc->reorder_block_Tensor( tmp, symm_orders[ii] );
           t_part->scale(symm_facs[ii]);
           cout << "t_part->norm() = " << t_part->norm() << endl;
           tensor_calc->put_sub_tensor( t_part, t_amps_full );
           
           break;
         }

         cout << "not found!!" << endl;
         if ( ii == 2 ) { 
           cout << "could not find symmetry block" << endl;
           throw logic_error( "could not build full t" ) ;
         }
       }
     }
   }

   return t_amps_full;     
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<SMITH::Tensor_<double>> CASPT2::CASPT2::build_full_H( shared_ptr<SMITH::Tensor_<double>> smith_h) {
cout << " CASPT2::CASPT2::build_full_H( shared_ptr<SMITH::Tensor_<double>>) " << endl;
/////////////////////////////////////////////////////////////////////////////

   auto tensor_calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>(); 

   shared_ptr<SMITH::Tensor_<double>> h_full = make_shared<SMITH::Tensor_<double>>(range_block_map_.at( "ffff") );
   h_full->allocate();
    
   vector<string> t_blocks = { "cvcv", "cvav", "avcv", "avav", "caca", "caaa", "aaca", "cacv", "caav", "aacv", "aaav", "cvca", "cvaa", "avca", "avaa" }; 

   vector<vector<int>> symm_orders = {{0,3,2,1},{2,1,0,3},{2,3,0,1}};  
   vector<double> symm_facs = { -1.0, -1,0, 1.0 };  
 
   auto reorder_vec = []( vector<SMITH::IndexRange> tbr, vector<int> ro ) { 
     vector<SMITH::IndexRange> nbr = { tbr[ro[0]] , tbr[ro[1]] , tbr[ro[2]], tbr[ro[3]] }; 
     return nbr; 
   }; 

   for ( auto block_name : t_blocks ) {
     cout<< "block_name = " << block_name << endl;
     vector<SMITH::IndexRange> t_block_ranges = range_block_map_.at(block_name);
     vector<SMITH::Index> t_block_index = { t_block_ranges[0].range(0), t_block_ranges[1].range(0), t_block_ranges[2].range(0), t_block_ranges[3].range(0) };

     if (smith_h->exists(t_block_index)) {
       shared_ptr<SMITH::Tensor_<double>> t_part = Tensor_Arithmetic_Utils::get_sub_tensor( smith_h, t_block_ranges );
       tensor_calc->put_sub_tensor( t_part, h_full );

     } else {
       cout << " block [ "; cout.flush(); for( auto idx : t_block_ranges) { cout << idx.size() << " "; cout.flush(); } cout << "] not in smith_h" << endl;

       for ( int ii = 0; ii != 3; ++ii ) {

         vector<SMITH::IndexRange> reord_t_block_ranges = reorder_vec( t_block_ranges, symm_orders[ii] ); 
         cout << "trying block [ "; cout.flush(); for( auto idx : t_block_ranges) { cout << idx.size() << " "; cout.flush(); } cout << "] not in smith_h" << endl;
         vector<SMITH::Index> reord_t_block_index = { reord_t_block_ranges[0].range(0), reord_t_block_ranges[1].range(0), reord_t_block_ranges[2].range(0), reord_t_block_ranges[3].range(0) };

         if ( smith_h->exists(reord_t_block_index) ) { 
           cout << "Found!!" << endl;
           shared_ptr<SMITH::Tensor_<double>> tmp = Tensor_Arithmetic_Utils::get_sub_tensor( smith_h, reord_t_block_ranges );
           auto t_part = tensor_calc->reorder_block_Tensor( tmp, symm_orders[ii] );
           t_part->scale(symm_facs[ii]);
           cout << "t_part->norm() = " << t_part->norm() << endl;
           tensor_calc->put_sub_tensor( t_part, h_full );
           
           break;
         }

         cout << "not found!!" << endl;
         if ( ii == 2 ) { 
           cout << "could not find symmetry block" << endl;
           throw logic_error( "could not build full t" ) ;
         }
       }
     }
   }


   return h_full;     
}

#endif
