#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_computer.h>
#include <src/util/f77.h>
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<Tensor_<double>>>>
Equation_Computer::Equation_Computer::get_gammas(int MM , int NN, string gamma_name){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  //Gets list of necessary gammas which are needed by truncating the gamma_range vector
  //e.g. [a,a,a,a,a,a] --> [a,a,a,a] --> [a,a]
  shared_ptr<vector<string>> gamma_ranges_str  = GammaMap->at(gamma_name)->id_ranges;

  auto gamma_ranges =  vector<shared_ptr<vector<IndexRange>>>(gamma_ranges_str->size()/2);
  for (int ii = 0 ; ii !=gamma_ranges.size(); ii++ ){ 
    shared_ptr<vector<string>> gamma_ranges_str_tmp = make_shared<vector<string>>(gamma_ranges_str->begin(), gamma_ranges_str->end()-ii*2);
    gamma_ranges[ii] = Get_Bagel_IndexRanges( gamma_ranges_str_tmp); 
  } 
  
  //Gamma data in vector format 
  shared_ptr<vector<shared_ptr<VectorB>>> gamma_data_vec = compute_gammas( MM, NN ) ;
  auto gamma_tensors = make_shared<vector<shared_ptr<Tensor_<double>>>>(0);
  for ( int ii = gamma_ranges.size()-1; ii != -1;  ii-- ) {
 
     shared_ptr<vector<int>> range_lengths  = make_shared<vector<int>>(0); 
     for (IndexRange idrng : *(gamma_ranges[ii]) )
       range_lengths->push_back(idrng.range().size()-1); 
     
     shared_ptr<Tensor_<double>> new_gamma_tensor= make_shared<Tensor_<double>>(*(gamma_ranges[ii]));
     new_gamma_tensor->allocate();
     
     shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(gamma_ranges[ii]->size(),0);  
     shared_ptr<vector<int>> mins = make_shared<vector<int>>(gamma_ranges[ii]->size(),0);  

     do {
       
       vector<Index> gamma_id_blocks(gamma_ranges[ii]->size());
       for( int jj = 0 ;  jj != gamma_id_blocks.size(); jj++)
         gamma_id_blocks[jj] =  gamma_ranges[ii]->at(jj).range(block_pos->at(jj));

       vector<int> range_sizes = get_sizes(gamma_id_blocks);
       shared_ptr<vector<int>> gamma_tens_strides = get_Tens_strides(range_sizes);  
       int gamma_block_size = accumulate( range_sizes.begin(), range_sizes.end(), 1, std::multiplies<int>() );
       int gamma_block_pos = inner_product( block_pos->begin(), block_pos->end(), gamma_tens_strides->begin(),  0); 

       unique_ptr<double[]> gamma_data_block(new double[gamma_block_size])  ;
       std::fill_n(gamma_data_block.get(), gamma_block_size, 0.0);
       blas::ax_plus_y_n(1.0,  gamma_data_vec->at(ii)->data()+gamma_block_pos, gamma_block_size, gamma_data_block.get());
       new_gamma_tensor->put_block( gamma_data_block, gamma_id_blocks);
     
     } while (fvec_cycle(block_pos, range_lengths, mins ));
     gamma_tensors->push_back(new_gamma_tensor);

     cout << "Printing Gamma of order " << gamma_ranges[ii]->size()/2  << endl;
     Print_Tensor(new_gamma_tensor);
  }
  cout << "out of loop" << endl;
   
  //Hack, fix this and use map instead
  reverse(gamma_tensors->begin(), gamma_tensors->end() ) ;
  
  return gamma_tensors;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma1  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<VectorB>>> 
Equation_Computer::Equation_Computer::compute_gammas(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12 MM = " << MM << " NN = " << NN  << endl;

  if (det_->compress()) { // uncompressing determinants
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }
  shared_ptr<Civec> ccbra = make_shared<Civec>(*cc_->data(MM));
  shared_ptr<Civec> ccket = make_shared<Civec>(*cc_->data(NN));

  shared_ptr<RDM<1>> gamma1x;
  shared_ptr<RDM<2>> gamma2x;
  shared_ptr<RDM<3>> gamma3x;
  tie(gamma1x, gamma2x, gamma3x) = compute_gamma12_from_civec(ccbra, ccket);
  
  auto  gamma1 = make_shared<VectorB>(norb_*norb_);
  copy_n(gamma1x->data(), norb_*norb_, gamma1->data());

  auto  gamma2 = make_shared<VectorB>(norb_*norb_*norb_*norb_);
  copy_n(gamma2x->data(), norb_*norb_*norb_*norb_, gamma2->data());

  auto  gamma3 = make_shared<VectorB>(norb_*norb_*norb_*norb_*norb_*norb_);
  copy_n(gamma3x->data(), norb_*norb_*norb_*norb_*norb_*norb_, gamma3->data());
  cc_->set_det(det_); 

  auto gamma_vec = make_shared<vector<shared_ptr<VectorB>>>(vector<shared_ptr<VectorB>> { gamma1, gamma2, gamma3}) ;
                    
 return gamma_vec;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the gamma matrix g_ij with elements c*_{M,I}< I | a*_{i} a_{j} | J > c_{N,J}
// mangled version of routines in fci_rdm.cc
// can use RDM type for convenience, but everything by gamma1  is _not_ an rdm 
///////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12(const int MM, const int NN ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12 MM = " << MM << " NN = " << NN  << endl;

  if (det_->compress()) { // uncompressing determinants
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }
  shared_ptr<Civec> ccbra = make_shared<Civec>(*cc_->data(MM));
  shared_ptr<Civec> ccket = make_shared<Civec>(*cc_->data(NN));
 
  shared_ptr<RDM<1>> gamma1;
  shared_ptr<RDM<2>> gamma2;
  shared_ptr<RDM<3>> gamma3;
  tie(gamma1, gamma2, gamma3) = compute_gamma12_from_civec(ccbra, ccket);
 
  cc_->set_det(det_); 

  return tie(gamma1, gamma2, gamma3);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "compute_gamma12_from_civec" << endl;

  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);
  
  shared_ptr<Dvec> dket;
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb_*norb_);
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }
  
  return compute_gamma12_last_step(dbra, dket, cbra);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>, shared_ptr<RDM<3>> >
Equation_Computer::Equation_Computer::compute_gamma12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "compute_gamma12_last_step" << endl;

  const int nri = cibra->asize()*cibra->lenb();
  const int n1 = norb_;  
  const int n2 = n1*norb_;  
 
  // gamma1 c^dagger <I|\hat{E}|0>
  // gamma2 \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto gamma1 = make_shared<RDM<1>>(norb_);
  auto gamma2 = make_shared<RDM<2>>(norb_);
  auto gamma3 = make_shared<RDM<3>>(norb_);

  //section to be made recursive for arbitrary orders of gamma
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, n2);
    for (int i = 0; i != n2; ++i)
      copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
 
    auto gamma1t = btas::group(*gamma1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, gamma1t, {1});

    auto dbra_data = dket_data;
    if (dbra != dket) {
      dbra_data = make_shared<Matrix>(nri, n2);
      for (int i = 0; i != n2; ++i)
        copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
    }

    const char   transa = 'N';
    const char   transb = 'T';
    const double alpha = 1.0;
    const double beta = 0.0; 
    const int n3 = n2*norb_;  
    const int n4 = n2*n2;
 
    //Very bad way of getting gamma3  [sum_{K}<I|i*j|K>.[sum_{L}<K|k*l|L>.<L|m*n|J>]]
    for ( int mn = 0; mn!=n2 ; mn++){
      Matrix gamma3_klmn(n2, n2);
      int m = mn/n1;
      int n = mn-(m*n1);
      auto dket_klmn = make_shared<Dvec>(dket->det(), n2);
      sigma_2a1(dket->data(mn), dket_klmn);
      sigma_2a2(dket->data(mn), dket_klmn);

      //make gamma3_klmn block
      auto dket_data_klmn = make_shared<Matrix>(nri, n2);
      for (int kl = 0; kl != n2; kl++)
        copy_n(dbra->data(kl)->data(), nri, dket_data_klmn->element_ptr(0, kl));
      btas::contract(1.0, *dbra_data, {1,0}, *dket_data_klmn, {1,2}, 0.0, gamma3_klmn, {0,2});
       
      //copy gamma3_klmn block into gamma3 and undo transpose of dbra in btas::contract
      unique_ptr<double[]> buf(new double[norb_*norb_]);
      for (int k = 0; k != n1; ++k) {
        for (int l = 0; l != n1; ++l) {
          copy_n(gamma3_klmn.data(), n2, buf.get());
          blas::transpose(buf.get(), n1, n1, gamma3->element_ptr(0,0,k,l,m,n));
        }
      }
    }

    auto gamma2t = group(group(*gamma2, 2,4), 0,2);
    btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, gamma2t, {0,2});
  }
 
  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&gamma2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, gamma2->element_ptr(0,0,k,i));
    }
  }
 
  return tie(gamma1, gamma2, gamma3);
}

//////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
//////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a1(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma) const {
//////////////////////////////////////////////////////////////////////////////////////////////
  //cout << "sigma_2a1" << endl;
  const int lb = sigma->lenb();
  const int ij = sigma->ij();
  const double* const source_base = cvec->data();

  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = sigma->data(ip)->data();

    for (auto& iter : cvec->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      blas::ax_plus_y_n(sign, source_base + iter.target*lb, lb, target_array);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Taken directly from src/ci/fci/knowles_compute.cc         
// I'm sure there was a version which used transposition of the civector; this looks slow.
///////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::sigma_2a2(shared_ptr<const Civec> cvec, shared_ptr<Dvec> sigma) const {
///////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "sigma_2a2" << endl;
  const int la = sigma->lena();
  const int ij = sigma->ij();

  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cvec->element_ptr(0, i);

    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = sigma->data(ip)->element_ptr(0, i);

      for (auto& iter : cvec->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_gamma_tensor( int MM , int NN, string gamma_name) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  if( Gamma_data_map->find(gamma_name) != Gamma_data_map->end()){ 

    cout << "already have data for " << gamma_name << endl;  
 
  } else { 

    //for now just use specialized routines, this must be made generic at some point
    if (GammaMap->at(gamma_name)->id_ranges->size() == 2 ) { 
      get_gamma_2idx( MM, NN, gamma_name ) ;
      cout << "------------------ "<<  gamma_name  << " ---------------------" << endl; 
      Print_Tensor(Gamma_data_map->at(gamma_name));
    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 4 ) { 
       //get rdm 4
    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 6 ) { 
       //get rdm 6
    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 8 ) { 
       //get rdm 8
    }    
  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_gamma_2idx(const int MM, const int NN, string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "get_gamma_2idx MM = " << MM << " NN = " << NN  << endl;
  cout << "gamma_name = " << gamma_name << endl;

  if (det_->compress()) { // uncompressing determinants, probably should be done in the sigma_blocked routine
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<Civec> ccbra = make_shared<Civec>(*cc_->data(MM));
  shared_ptr<Civec> ccket = make_shared<Civec>(*cc_->data(NN));

  // build gamma tensor
  shared_ptr<vector<IndexRange>> gamma_ranges  = Get_Bagel_IndexRanges( GammaMap->at(gamma_name)->id_ranges );   
  shared_ptr<vector<int>>        range_lengths = get_range_lengths(gamma_ranges); 

  shared_ptr<Tensor_<double>>    gamma_tensor  = make_shared<Tensor_<double>>(*(gamma_ranges));
  gamma_tensor->allocate();
  gamma_tensor->zero();
  
  shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(gamma_ranges->size(),0);  
  shared_ptr<vector<int>> mins = make_shared<vector<int>>(gamma_ranges->size(),0);  
  
  shared_ptr<vector<pair<size_t,size_t>>> gamma_block_start;

  do {
    
    vector<Index> gamma_id_blocks = *(get_rng_blocks( block_pos, *gamma_ranges));
  
    shared_ptr<vector<pair<size_t,size_t>>> ij_interval = get_block_start( gamma_ranges, block_pos ) ;

    unique_ptr<double[]> gamma_data_block = gamma_2idx_block( ccbra, ccket, ij_interval->at(0), ij_interval->at(1) );

    gamma_tensor->put_block( gamma_data_block, gamma_id_blocks);
  
  } while (fvec_cycle(block_pos, range_lengths, mins ));

  Gamma_data_map->emplace(gamma_name, gamma_tensor); 

  cc_->set_det(det_); 

  return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unique_ptr<double[]>
Equation_Computer::Equation_Computer::gamma_2idx_block( shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket,
                                                        pair<size_t,size_t> irange, pair<size_t,size_t> jrange     ) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "gamma_2idx_block" << endl;

  unique_ptr<double[]> sigma_block = sigma_blocked( cbra, irange, jrange );
  size_t iblock_size  = irange.second - irange.first; 
  size_t jblock_size  = jrange.second - jrange.first; 

  unique_ptr<double[]> gamma_block(new double[iblock_size*jblock_size])  ;
  std::fill_n(gamma_block.get(), iblock_size*jblock_size, 0.0);
  
  dgemv_("N", cbra->size(), iblock_size*jblock_size, 1.0, sigma_block.get(), cbra->size(), cket->data(), 1, 0.0, gamma_block.get(), 1);
  
  return gamma_block;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unique_ptr<double[]> 
Equation_Computer::Equation_Computer::sigma_blocked(shared_ptr<const Civec> cvec, 
                                                    pair<size_t,size_t> irange, pair<size_t,size_t> jrange) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "sigma_blocked" << endl;

  const size_t lena = cvec->lena();
  const size_t lenb = cvec->lenb();

  size_t iblock_size  = irange.second - irange.first; 
  size_t jblock_size  = jrange.second - jrange.first; 

  const size_t sigma_block_size =lena*lenb*(iblock_size*norb_+jblock_size);

  unique_ptr<double[]> sigma_block(new double[sigma_block_size])  ;
  std::fill_n(sigma_block.get(), sigma_block_size, 0.0);

  const double* const source_base = cvec->data();
  for (size_t ii = 0; ii != iblock_size; ++ii) {
    for (size_t jj = 0; jj != jblock_size; ++jj) {
      double* const target_base = sigma_block.get() + ((ii*norb_+jj)*lena*lenb);
      
      for (auto& iter : cvec->det()->phia( (ii+irange.first)*norb_ + (jj+jrange.first) )) {
        const double sign = static_cast<double>(iter.sign);
        double* const target_array = target_base + iter.source*lenb;
        blas::ax_plus_y_n(sign, source_base + iter.target*lenb, lenb, target_array);
      }
    }
  }

  //should do this by transposing civector so can use blas copy as above
  for (int adet_pos = 0; adet_pos < lena; ++adet_pos) {
    const double* const source_array0 = cvec->element_ptr(0, adet_pos);

    for (size_t ii = 0; ii != iblock_size; ++ii) {
      for (size_t jj = 0; jj != jblock_size; ++jj) {
        double* const target_array0 = sigma_block.get() + ((ii*norb_+jj)*lena*lenb) + adet_pos;

        for (auto& iter : cvec->det()->phib( (ii+irange.first)*norb_+(jj+jrange.first) ) ) {
          const double sign = static_cast<double>(iter.sign);
          target_array0[iter.source] += sign * source_array0[iter.target];
        }
      }
    }    
  }      
  
  return sigma_block;

}        

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Out puts the interval [start_orb, end_orb) of the index blocks at block_pos
// Note the interval is relative to the input IndexRanges.
////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<size_t, size_t>>>
Equation_Computer::Equation_Computer::get_block_start( shared_ptr<vector<IndexRange>> id_ranges, 
                                                       shared_ptr<vector<int>> block_pos         ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "Equation_Computer::get_block_info" << endl;

  vector<pair<size_t,size_t>> block_start_end(block_pos->size());
  for (int ii = 0 ; ii != block_start_end.size() ; ii++){
    size_t block_start = 0;

    for (int jj = 0 ; jj != block_pos->at(ii); jj++) 
      block_start += id_ranges->at(ii).range(jj).size();

    block_start_end[ii] = make_pair(block_start, block_start+id_ranges->at(ii).range(block_pos->at(ii)).size());

  }
//  cout << "block_start_end = [ " ; cout.flush(); for ( auto elem : block_start_end  ) {cout << "(" << elem.first <<  ","  << elem.second << ") " ; } cout << " ] " << endl;

  return make_shared<vector<pair<size_t,size_t>>>(block_start_end);
}

#endif
