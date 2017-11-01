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
Equation_Computer::Equation_Computer::get_gammas(int MM, int NN, string gamma_name){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  //Gets list of necessary gammas which are needed by truncating the gamma_range vector
  //e.g. [a,a,a,a,a,a] --> [a,a,a,a] --> [a,a]
  shared_ptr<vector<string>> gamma_ranges_str = GammaMap->at(gamma_name)->id_ranges;

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

  string det_name = get_civec_name( MM, cc_->data(MM)->det()->norb(),  cc_->data(MM)->det()->nelea(), cc_->data(MM)->det()->neleb());  
  shared_ptr<const Determinants> det_ =  determinants_map->at(det_name);
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

  shared_ptr<const Determinants> det_ =  determinants_map->at(get_det_name(cc_->data(MM)->det()));
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

    //DetMap(size_t t, int si, size_t s, int o) : target(t), sign(si), source(s), ij(o) {}
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
      build_gamma_2idx_tensor( MM, NN, nelea_, neleb_, norb_,  gamma_name ) ;
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

  string det_name = get_civec_name( MM, cc_->data(MM)->det()->norb(),  cc_->data(MM)->det()->nelea(), cc_->data(MM)->det()->neleb());  
  shared_ptr<const Determinants> det_ =  determinants_map->at(det_name);
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
  
  do {
    
    vector<Index> gamma_id_blocks = *(get_rng_blocks( block_pos, *gamma_ranges));
  
    shared_ptr<vector<pair<size_t,size_t>>> ij_interval = get_block_start( gamma_ranges, block_pos ) ;

    unique_ptr<double[]> gamma_data_block = gamma_2idx_block( ccbra, ccket, ij_interval->at(0), ij_interval->at(1) );

    gamma_tensor->put_block( gamma_data_block, gamma_id_blocks);
  
  } while (fvec_cycle(block_pos, range_lengths, mins ));

  Gamma_data_map->emplace(gamma_name, gamma_tensor); 
  Data_map->emplace(gamma_name, gamma_tensor); 

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
  cout << "sigma_block_size = " <<  sigma_block_size << endl;
  cout << "irange.second , irange.first = "<<  irange.second << " , " <<  irange.first  << endl;
  cout << "jrange.second , jrange.first = "<<  jrange.second << " , " <<   jrange.first  << endl;
  cout << "lena , lenb                  = "<<  lena  << " , " <<  lenb  << endl;

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> 
Equation_Computer::Equation_Computer::convert_civec_to_tensor( shared_ptr<const Civec> civector, int state_num ) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::convert_civec_to_tensor" << endl;

  //NOTE: must be adapted to handle arbitrary spin sectors
  string civec_name = get_civec_name(state_num, civector->det()->norb(), civector->det()->nelea(), civector->det()->neleb());  
  vector<IndexRange> civec_idxrng(1, *(range_conversion_map->at(civec_name)) );  

  cout <<" civec_idxrng[0].nblock()       = " << civec_idxrng[0].nblock()     <<  endl;
  cout <<" civec_idxrng[0].size()         = " << civec_idxrng[0].size()       <<  endl;
  cout <<" civec_idxrng[0].range().size() = " << civec_idxrng[0].range().size() <<  endl;
  cout <<" civec_idxrng[0].range().size() = " << civec_idxrng[0].range().size() <<  endl;
  
  shared_ptr<Tensor_<double>> civec_tensor = make_shared<Tensor_<double>>( civec_idxrng );
  civec_tensor->allocate();
  civec_tensor->zero();

  size_t idx_position = 0;

  cout << "civectordata = " ; cout.flush(); 
  for ( Index idx_block : civec_idxrng[0].range() ){
     unique_ptr<double[]> civec_block(new double[idx_block.size()]);
     std::fill_n(civec_block.get(), idx_block.size(), 0.0);
     copy_n( civector->data() + idx_position, idx_block.size(), civec_block.get());

     for ( int ii = 0 ; ii != idx_block.size() ; ii++ ) 
       cout << *(civector->data() + idx_position + ii) << " "; 
     cout.flush();
  
     civec_tensor->add_block(civec_block, vector<Index>({ idx_block })) ;  
     idx_position += idx_block.size();  
  }

  cout <<endl;

  //will have to modify for relativistic case
  CIvec_data_map->emplace( civec_name, civec_tensor); 
  determinants_map->emplace( civec_name, civector->det() ); 

  Data_map->emplace( civec_name, civec_tensor); 

  return civec_tensor;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
string Equation_Computer::Equation_Computer::get_civec_name(int state_num, int norb, int nalpha, int nbeta) const { 
////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = to_string(state_num) + "_["+ to_string(norb)+"{" + to_string(nalpha) + "a," + to_string(nbeta) + "b}]" ;
  return name ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Outputs the interval [start_orb, end_orb) of the index blocks at block_pos
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Equation_Computer::Equation_Computer::build_gamma_2idx_tensor( int NN, int MM, int nelea, int neleb, int norb, string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "build_gamma_2idx_tensor : " << gamma_name << endl;

  shared_ptr<vector<string>> gamma_ranges_str  = GammaMap->at(gamma_name)->id_ranges;

  string Bra_name = get_civec_name( NN, norb, nelea, neleb );  
  string Ket_name = get_civec_name( MM, norb, nelea, neleb );  

  shared_ptr<vector<bool>> aops = make_shared<vector<bool>>(vector<bool> { true, false } );
  string sigma_name = get_sigma_name( Bra_name, Ket_name, gamma_ranges_str, aops );

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
     
     cout << "already got " << sigma_name << ", so use it " << endl;
  
  } else { //should replace this with blockwise and immediate contraction build.

    build_sigma_2idx_tensor( Bra_name, Ket_name, gamma_ranges_str);
   
    shared_ptr<Tensor_<double>> gamma_2idx = contract_different_tensors( make_pair(0,0), Bra_name, sigma_name, gamma_name);
     
    Gamma_data_map->emplace( gamma_name, gamma_2idx );
    Data_map->emplace( gamma_name, gamma_2idx );
 
    Print_Tensor(Gamma_data_map->at(gamma_name));

  }
 
  return;
 
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Equation_Computer::Equation_Computer::build_sigma_2idx_tensor(string Bra_name, string Ket_name, shared_ptr<vector<string>> orb_ranges_str)  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_sigma_2idx_tensor" << endl;
 
  // temp hack for aops, must implement SigmaInfo class and map so this is cleaner
  shared_ptr<vector<bool>> aops = make_shared<vector<bool>>(vector<bool> { true, false } );
  string sigma_name = get_sigma_name( Bra_name, Ket_name, orb_ranges_str, aops );

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

  } else { 
    
    shared_ptr<Tensor_<double>> Bra_civec = CIvec_data_map->at(Bra_name);
    shared_ptr<const Determinants> Bra_det = determinants_map->at(Bra_name);
    
    shared_ptr<vector<IndexRange>> orb_ranges  = Get_Bagel_IndexRanges( orb_ranges_str );   
    
    shared_ptr<vector<IndexRange>> sigma_ranges = make_shared<vector<IndexRange>>( Bra_civec->indexrange() );
    sigma_ranges->insert(sigma_ranges->end(), orb_ranges->begin(), orb_ranges->end());
    
    shared_ptr<Tensor_<double>> sigma_tensor = make_shared<Tensor_<double>>(*(sigma_ranges));
    sigma_tensor->allocate();
    sigma_tensor->zero();
   
    //note that we because of the loop ranges has 
    shared_ptr<vector<int>> range_lengths = get_range_lengths(sigma_ranges); 
    shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(sigma_ranges->size(),0);  
    shared_ptr<vector<int>> mins = make_shared<vector<int>>(sigma_ranges->size(),0);  
    vector<int> sigma_offsets(3,0);
    
    // loop through sigma ranges,  loop over Ket ranges inside build_sigma_block;  sigma is much larger
    // than Ket, so don't want to move it about.
    do {
       
      vector<Index> sigma_id_blocks = *(get_rng_blocks( block_pos, *sigma_ranges));
      build_sigma_block( sigma_name, sigma_id_blocks, sigma_offsets, Ket_name ) ;
      
      for (int ii = 0 ; ii != 3; ii++ )
        sigma_offsets[ii] += sigma_id_blocks[ii].size();
    
    } while (fvec_cycle(block_pos, range_lengths, mins ));
    
    Sigma_data_map->emplace(sigma_name, sigma_tensor); 
    Data_map->emplace(sigma_name, sigma_tensor); 

  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::build_sigma_block( string sigma_name, vector<Index>& sigma_id_blocks,
                                                              vector<int>& sigma_offsets, string Ket_name  ) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::build_sigma_block" << endl; 

  size_t Bra_block_size = sigma_id_blocks[0].size();
  size_t iblock_size    = sigma_id_blocks[1].size();
  size_t jblock_size    = sigma_id_blocks[2].size();

  const size_t sigma_block_size = Bra_block_size*(iblock_size*norb_+jblock_size);

  cout << " iblock_size      = " << iblock_size << endl;
  cout << " jblock_size      = " << jblock_size << endl;
  cout << " Bra_block_size   = " << Bra_block_size << endl;
  cout << " sigma_block_size = " << sigma_block_size << endl;
  unique_ptr<double[]> sigma_block(new double[sigma_block_size])  ;
  std::fill_n(sigma_block.get(), sigma_block_size, 0.0);

  shared_ptr<const Determinants> Ket_det = determinants_map->at(Ket_name);
  shared_ptr<Tensor_<double>> Ket_Tens  = CIvec_data_map->at(Ket_name);
  shared_ptr<IndexRange> Ket_range       = range_conversion_map->at(Ket_name);
  // size_t lb = Ket_det->lenb() <= Ket_idx_block.size() ? Ket_det_lenb() : Ket_idx_block.size();      
  size_t Ket_offset = 0;

  auto in_Bra_range = [&sigma_offsets, &Bra_block_size ]( size_t& pos) {
      if (( pos > (sigma_offsets[0]+Bra_block_size)) || ( pos < sigma_offsets[0] ) ){
         cout << "out of sigma range " << endl;
         cout << "pos + len      = " << pos << endl;
         cout << "Bra_offset     = " << sigma_offsets[0] << endl;
         cout << "Bra_Block_size = " << Bra_block_size << endl;
         return false;
      }
      return true;
      };

   auto in_Ket_range = [&Ket_offset ]( size_t& pos , size_t Ket_block_size ) {
      if ( ( pos  > (Ket_offset+ Ket_block_size )) || ( pos < Ket_offset ) ){
         cout << "out of Ket range " << endl;
         cout << "pos             = " << pos << endl;
         cout << "Ket_offset      = " << Ket_offset << endl;
         cout << "Ket_Block_size  = " << Ket_block_size << endl;
         return false;
      }
      return true;
      };

   const int lena = Ket_det->lena();
   const int lenb = Ket_det->lenb();
  
  // vector<DetMap> phia ,  DetMap(size_t t, int si, size_t s, int o) : target(t), sign(si), source(s), ij(o) {}
  // source and target seem backward, but it's this way round everywhere else, strange name perhaps?
  // Must be fixed to use BLAS, however, you will have to be careful about the ranges; 
  // ci_block needs to have a length which is some integer multiple of lenb.
 
  double* sigma_ptr = sigma_block.get();
  for ( Index Ket_idx_block : Ket_range->range()) { 
    unique_ptr<double[]> Ket_block = Ket_Tens->get_block(vector<Index>{Ket_idx_block});
    double* Ket_ptr = Ket_block.get();
    size_t Ket_block_size = Ket_idx_block.size();
    
    for( size_t ii = sigma_offsets[1]; ii != iblock_size + sigma_offsets[1]; ++ii) {
      for( size_t jj = sigma_offsets[2]; jj != iblock_size + sigma_offsets[2]; ++jj) {
        for (auto& iter : Ket_det->phia(ii, jj)) {
          size_t Sshift = iter.source*lenb;
          size_t Kshift = iter.source*lenb;
          double* sigma_pos = sigma_block.get() + iter.source *lenb;
          double* Ket_pos   = Ket_block.get() +   iter.target *lenb;
          double sign = static_cast<double>(iter.sign);
          for( size_t ib = 0; ib != lenb; ++ib, sigma_pos++, Ket_pos++, Sshift++, Kshift++) {
            if( (in_Ket_range(Kshift, Ket_block_size) && in_Bra_range(Sshift) ))
              *sigma_pos += (*Ket_pos * sign);
          }
        }
      }
    }

    Ket_ptr = Ket_block.get();
    double* sigma_block_ia = sigma_block.get();
    for (int ia = 0; ia < lena; ++ia) {
      Ket_ptr += lenb;
      sigma_block_ia += lenb;
      for( size_t ii = sigma_offsets[1]; ii != iblock_size + sigma_offsets[1]; ++ii) {
        for( size_t jj = sigma_offsets[2]; jj != iblock_size + sigma_offsets[2]; ++jj) {
          for (auto& iter : Ket_det->phib(ii,jj)) {
            const double sign = static_cast<double>(iter.sign);
            size_t shift1 = ia*lenb+iter.target;
            size_t shift2 = ia*lenb+iter.source; 
            if ( (in_Ket_range( shift1,  Ket_block_size)) && (in_Bra_range( shift2 )) )
              *(sigma_block.get() + iter.source) +=  *(Ket_ptr + iter.target) * sign ;
          }
        }
      }
    }
    Ket_offset += Ket_idx_block.size(); 
  }
  
  return;
 
}       
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string Equation_Computer::Equation_Computer::get_sigma_name( string Bra_name, string Ket_name , shared_ptr<vector<string>>  orb_ranges,
                                                             shared_ptr<vector<bool>>  aops ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   string sigma_name = Bra_name + "_";

   for (string rng : *orb_ranges ) 
     sigma_name += rng[0];

   for ( bool aop : *aops ) {
     if ( aop ) {
       sigma_name += "1";
     } else { 
       sigma_name += "0";
     }
   } 

  sigma_name += "_"+Ket_name;

  return sigma_name;
}
 
#endif
