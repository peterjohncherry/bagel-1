#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/gamma_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Gamma_Computer::Gamma_Computer::Gamma_Computer( shared_ptr<Equation<double>> eqn_info_in,
						shared_ptr<map< string, shared_ptr<IndexRange>>> range_conversion_map_in ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::Gamma_Computer" << endl;
  maxtile = 10000;
  GammaMap = eqn_info_in->GammaMap;
  
  CIvec_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();
  Sigma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();
  Gamma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();
  Determinants_map = make_shared<std::map<std::string, std::shared_ptr<const Determinants>>>();

  cimaxblock = 100; //figure out what is best, maxtile is 10000, so this is chosen to have one index block. Must be consistent if contraction routines are to work...

  Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>();

  tester();
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::Initialize_wfn_info( shared_ptr<Civec> civector, int state_num ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  nelea_ = civector->det()->nelea();
  neleb_ = civector->det()->neleb();
  norb_  = civector->det()->norb();
  get_civector_indexranges(1);

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::get_gamma_tensor( int MM , int NN, string gamma_name) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::get_gammas"  << endl;
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
void Gamma_Computer::Gamma_Computer::get_gamma_2idx(const int MM, const int NN, string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "get_gamma_2idx MM = " << MM << " NN = " << NN  << endl;
  cout << "gamma_name = " << gamma_name << endl;

  string det_name = get_civec_name( MM, cc_->data(MM)->det()->norb(),  cc_->data(MM)->det()->nelea(), cc_->data(MM)->det()->neleb());  
  shared_ptr<const Determinants> det_ =  Determinants_map->at(det_name);
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

  cc_->set_det(det_); 

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
unique_ptr<double[]>
Gamma_Computer::Gamma_Computer::gamma_2idx_block( shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket,
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
Gamma_Computer::Gamma_Computer::sigma_blocked( shared_ptr<const Civec> cvec,
                                               pair<size_t,size_t> irange, pair<size_t,size_t> jrange ) const {
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
Gamma_Computer::Gamma_Computer::convert_civec_to_tensor( shared_ptr<const Civec> civector, int state_num ) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::convert_civec_to_tensor" << endl;

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
  Determinants_map->emplace( civec_name, civector->det() ); 


  return civec_tensor;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_gamma_2idx_tensor( int NN, int MM, int nelea, int neleb, int norb, string gamma_name ) {
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
   
    shared_ptr<Tensor_<double>> gamma_2idx = Tensor_Calc->contract_different_tensors(  CIvec_data_map->at(Bra_name), Sigma_data_map->at(sigma_name),  make_pair(0,0));
     
    Gamma_data_map->emplace( gamma_name, gamma_2idx );
 
    Print_Tensor(Gamma_data_map->at(gamma_name));

  }
 
  return;
 
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_sigma_2idx_tensor(string Bra_name, string Ket_name, shared_ptr<vector<string>> orb_ranges_str)  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_sigma_2idx_tensor" << endl;
 
  // temp hack for aops, must implement SigmaInfo class and map so this is cleaner
  shared_ptr<vector<bool>> aops = make_shared<vector<bool>>(vector<bool> { true, false } );
  string sigma_name = get_sigma_name( Bra_name, Ket_name, orb_ranges_str, aops );

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

  } else { 
    
    shared_ptr<Tensor_<double>> Bra_civec = CIvec_data_map->at(Bra_name);
    shared_ptr<const Determinants> Bra_det = Determinants_map->at(Bra_name);
    
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

  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigma_block( string sigma_name, vector<Index>& sigma_id_blocks,
                                                        vector<int>& sigma_offsets, string Ket_name  ) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigma_block" << endl; 

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

  shared_ptr<const Determinants> Ket_det = Determinants_map->at(Ket_name);
  shared_ptr<Tensor_<double>> Ket_Tens  = CIvec_data_map->at(Ket_name);
  shared_ptr<IndexRange> Ket_range       = range_conversion_map->at(Ket_name);
  // size_t lb = Ket_det->lenb() <= Ket_idx_block.size() ? Ket_det_lenb() : Ket_idx_block.size();      
  size_t Ket_offset = 0;

  auto in_Bra_range = [&sigma_offsets, &Bra_block_size ]( size_t& pos) {
      if (( pos > (sigma_offsets[0]+Bra_block_size)) || ( pos < sigma_offsets[0] ) ){
         cout << "out of sigma range " << endl;
         cout << "pos            = " << pos << endl;
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
  
  // Must be changed to use BLAS, however, you will have to be careful about the ranges; 
  // ci_block needs to have a length which is some integer multiple of lenb.
  double* sigma_ptr = sigma_block.get();
  for ( Index Ket_idx_block : Ket_range->range()) { 
    unique_ptr<double[]> Ket_block = Ket_Tens->get_block(vector<Index>{Ket_idx_block});
    double* Ket_ptr = Ket_block.get();
    size_t Ket_block_size = Ket_idx_block.size();
    
    for( size_t ii = sigma_offsets[1]; ii != iblock_size + sigma_offsets[1]; ++ii) {
      for( size_t jj = sigma_offsets[2]; jj != iblock_size + sigma_offsets[2]; ++jj) {
        for ( DetMap iter : Ket_det->phia(ii, jj)) {
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
          for ( DetMap iter : Ket_det->phib(ii,jj)) {
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//should be extended to deal with spin sectors
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::get_civector_indexranges(int nstates) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  for ( int ii = 0 ; ii != nstates; ii++ ) 
    range_conversion_map->emplace( get_civec_name( ii , cc_->data(ii)->det()->norb(), cc_->data(ii)->det()->nelea(), cc_->data(ii)->det()->neleb()),
                                   make_shared<IndexRange>(cc_->data(ii)->lena()*cc_->data(ii)->lenb(), cimaxblock ));  

  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string Gamma_Computer::Gamma_Computer::get_sigma_name( string Bra_name, string Ket_name , shared_ptr<vector<string>>  orb_ranges,
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
////////////////////////////////////////////////////////////////////////////////////////////////////////
string Gamma_Computer::Gamma_Computer::get_civec_name(int state_num, int norb, int nalpha, int nbeta) const { 
////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = to_string(state_num) + "_["+ to_string(norb)+"{" + to_string(nalpha) + "a," + to_string(nbeta) + "b}]" ;
  return name ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string Gamma_Computer::Gamma_Computer::get_det_name(shared_ptr<const Determinants> Detspace ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 return "[" +to_string(Detspace->norb()) + ",{"+to_string(Detspace->nelea())+"a,"+to_string(Detspace->neleb())+"b}]";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::tester(){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
  shared_ptr<Tensor_<double>> civec1 =  convert_civec_to_tensor( cc_->data(0), 0 );
  shared_ptr<Tensor_<double>> civec2 =  convert_civec_to_tensor( cc_->data(0), 0 );

  double normval = civec1->dot_product(civec2); 
  cout << " civec1->dot_product(civec2) = " << normval << endl;
  cout << " civec1->rms()               = " << civec1->rms()  << endl;
  cout << " civec1->norm()              = " << civec1->norm() << endl;
  
  assert(!(abs(normval -1.00) > 0.000000000001) ); 
  
  return;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> Gamma_Computer::Gamma_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Gamma_Computer::Get_Bagel_IndexRanges 1arg" << endl;

  shared_ptr<vector<IndexRange>> ranges_Bagel = make_shared<vector<IndexRange>>(ranges_str->size());
  for ( int ii =0 ; ii != ranges_str->size(); ii++) 
    ranges_Bagel->at(ii) = *range_conversion_map->at(ranges_str->at(ii));

  return ranges_Bagel;
}

#endif
