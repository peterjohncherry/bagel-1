#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/b_gamma_computer.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
B_Gamma_Computer::B_Gamma_Computer::B_Gamma_Computer( std::shared_ptr<const Dvec > cc_in, 
                                                      std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                                                      std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> Gamma_info_map_in,
                                                      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Gamma_data_map_in ,
                                                      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Sigma_data_map_in ,
                                                      std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CIvec_data_map_in ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::B_Gamma_Computer::B_Gamma_Computer" << endl;
  cc_     = cc_in;
  
  Gamma_info_map = Gamma_info_map_in;
  Gamma_data_map = Gamma_data_map_in;
  Sigma_data_map = Sigma_data_map_in;
  CIvec_data_map = CIvec_data_map_in;

  range_conversion_map = range_conversion_map_in;

  Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>();

  dvec_sigma_map = make_shared<std::map< std::string, std::shared_ptr<Dvec>>>();
  det_old_map    = make_shared<std::map< std::string, std::shared_ptr<Determinants>>>();
  cvec_old_map   = make_shared<std::map< std::string, std::shared_ptr<Civec>>>();
}  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::get_gamma( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::get_gamma :  " <<  gamma_name << endl;

  if ( Gamma_data_map->find( gamma_name ) == Gamma_data_map->end()) {

    shared_ptr<GammaInfo> gamma_info = Gamma_info_map->at(gamma_name);

    //note this has reverse iterators!
    if (gamma_info->order > 2 ) { 
    
      get_gamma(gamma_info->prev_gamma_name());
      compute_sigmaN( gamma_info );
      convert_Dvec_sigma_to_tensor( gamma_info );
      get_gammaN_from_sigmaN ( gamma_info ) ;
    
    } else {
    
      compute_sigma2( gamma_info );
      convert_Dvec_sigma_to_tensor( gamma_info );
      get_gamma2_from_sigma2( gamma_info );
    
    }

  }

  return;                              
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Obviously a hack for now; duplicates everything.
// Assumes that cimaxblock is _larger_ than the length of the civec. 
// With cimaxblock = 10,000 this is OK up to (8,8), so big enough for lanthanides
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::convert_Dvec_sigma_to_tensor( shared_ptr<GammaInfo> gamma_info ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::Convert_to_sigma_dvec_tensor :" << gamma_info->sigma_name <<  endl;

  string sigma_name = gamma_info->sigma_name;
  int order = gamma_info->order;
  shared_ptr<Dvec> Dvec_sigma = dvec_sigma_map->at( sigma_name ) ;

  shared_ptr<vector<IndexRange>> sigma_ranges = Get_Bagel_IndexRanges( gamma_info->sigma_id_ranges ); 
  print_vector( *(gamma_info->sigma_id_ranges) , "sigma_id_ranges") ;  cout << endl;
  shared_ptr<Tensor_<double>> Tens_sigma = make_shared<Tensor_<double>>( *sigma_ranges ); 
  Tens_sigma->allocate();
  Tens_sigma->zero();
  Sigma_data_map->emplace( sigma_name, Tens_sigma ); 

  shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( *sigma_ranges ) ;

  shared_ptr<vector<int>> range_lengths  =  get_range_lengths( sigma_ranges ) ;
  shared_ptr<vector<int>> block_pos      =  make_shared<vector<int>>(range_lengths->size(),0);  
  shared_ptr<vector<int>> mins           =  make_shared<vector<int>>(range_lengths->size(),0);  

  vector<int> orb_id_strides(order);
  for ( int ii = 0 ; ii != order ; ii++ ) 
    orb_id_strides[order - ii - 1] = (int)pow(Dvec_sigma->det()->norb(), ii);

  auto boffset = [&block_offsets, &block_pos]( int pos ){ return block_offsets->at(pos).at(block_pos->at(pos)); }; 
  do {
      
    print_vector( *block_pos, "block_pos" ) ; cout<< endl; cout <<" order = " <<  order << endl;

    vector<Index> sigma_id_blocks = *(get_rng_blocks( block_pos, *sigma_ranges ));
    int block_size = 1 ;
    int orb_pos = 0;
    for ( int  ii = 0 ; ii != order+1; ii++ ){
      orb_pos += boffset(ii)*orb_id_strides[ii]; 
      block_size *= sigma_id_blocks[ii].size(); 
    }    

    unique_ptr<double[]> sigma_block( new double[block_size] );
    double* sigma_block_ptr = sigma_block.get();

    int ci_block_size = sigma_id_blocks[order].size();
    cout << " block_size = " << block_size << endl; cout << " ci_block_size = " << ci_block_size << endl;   cout << " block_size/ci_block_size = " << block_size/ci_block_size << endl; 
    for (int ii = 0; ii != block_size/ci_block_size ; ii++ ) { 
      std::copy( Dvec_sigma->data(orb_pos)->data()+boffset(order), Dvec_sigma->data(orb_pos)->data()+boffset(order)+ci_block_size, sigma_block_ptr );   
      sigma_block_ptr+=ci_block_size;
      orb_pos++;
    }

    Tens_sigma->put_block( sigma_block, sigma_id_blocks ) ; 

  if (order == 4) { 
  int n1 = Dvec_sigma->det()->norb();
  int n2 = n1*n1;
  int n3 = n1*n2;
  int cisize = Dvec_sigma->det()->size();

 
  cout << endl << endl;
  cout << " Sigma4 test dvec " << endl;
  for ( int xx =0 ; xx!= cisize; xx++ ) {
    cout << xx << " : ";
    for ( int qq = 0 ; qq != n1; qq++){
      for ( int rr = 0 ; rr != n1; rr++){
        for ( int ss = 0 ; ss != n1; ss++){
          for ( int tt = 0 ; tt != n1; tt++){
            cout << *( Dvec_sigma->data(qq*n3 + rr*n2 + ss*n1 + tt )->data()+xx) << " " ;  
          }
        }
      }
    }
    cout << endl;
  }
  cout << endl << endl;
  cout << "sigma4_ in Tensor column format" << endl;
  for ( int xx =0 ; xx!= cisize; xx++ ) {
    cout << xx << " : ";
    for ( int qq = 0 ; qq != n1; qq++){
      for ( int rr = 0 ; rr != n1; rr++){
        for ( int ss = 0 ; ss != n1; ss++){
          for ( int tt = 0 ; tt != n1; tt++){
            cout << sigma_block[ (qq*n3 + rr*n2 + ss*n1 + tt)*cisize + xx ] << " " ;
          }
        }
      }
    }
    cout << endl;
  }
  cout << endl << endl;
  }
  } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::get_gammaN_from_sigmaN( shared_ptr<GammaInfo> gammaN_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::get_gammaN_from_sigmaN " << endl;

  shared_ptr<Dvec> sigmaN = dvec_sigma_map->at(gammaN_info->sigma_name);
  string Bra_name = gammaN_info->Bra_name();
  shared_ptr<Civec> Bra = cvec_old_map->at( Bra_name );      
  shared_ptr<Determinants> Bra_det = det_old_map->at( Bra_name ) ; 

  int orb_dim = pow(Bra_det->norb(), (gammaN_info->order - 2) );
  int orb2    = Bra_det->norb()*Bra_det->norb();

  unique_ptr<double[]> gammaN( new double[orb_dim*orb2]);
  std::fill_n(gammaN.get(), orb_dim*orb2, 0.0);

  for ( int  ii = 0; ii != orb_dim*orb2; ii++) 
    gammaN[ii] = ddot_( Bra_det->size(), sigmaN->data(ii)->data(), 1, Bra->data(), 1); 

  shared_ptr<Tensor_<double>> Tens_sigmaN = Sigma_data_map->at(gammaN_info->sigma_name);
  convert_civec_to_tensor( Bra_name );
  shared_ptr<vector<IndexRange>> gamma_ranges = Get_Bagel_IndexRanges( gammaN_info->id_ranges ); 
  shared_ptr<Tensor_<double>> Tens_gammaN = Tensor_Calc->contract_tensor_with_vector( Tens_sigmaN, CIvec_data_map->at(Bra_name), make_pair(gammaN_info->order, 0) );
  Gamma_data_map->emplace( gammaN_info->name, Tens_gammaN ); 
 
  cout << endl; 
  Print_Tensor( Tens_gammaN , gammaN_info->name ) ; cout << endl << endl;


  if ( gammaN_info->order == 4 ) { 
    int norb  = Bra_det->norb();
    cout << "printing gamma4 from old method " << endl;
    int pos = 0 ; 
    for ( int ii = 0 ; ii != norb ; ii++) {
      for ( int jj = 0 ; jj != norb ; jj++) { 
        cout << " [ " << ii << " " << jj << " 0 0 ] "<< endl;
        for ( int kk = 0 ; kk != norb ; kk++) { 
          for ( int ll = 0 ; ll != norb ; ll++) { 
             cout << gammaN[pos++] << " " ;
          }
          cout << endl; 
        }
      }
    }

    shared_ptr<Dvec> sigma2 = dvec_sigma_map->at(gammaN_info->prev_sigma_name());
    unique_ptr<double[]> gamma2( new double[orb2]);
    for ( int  ii = 0; ii != orb2; ii++) 
      gamma2[ii] = ddot_( Bra_det->size(), sigma2->data(ii)->data(), 1, Bra->data(), 1); 

    cout << " TEST rdm4 " << endl;
    pos = 0;
    for ( int ii = 0 ; ii != norb ; ii++) {
      for ( int jj = 0 ; jj != norb ; jj++) { 
        cout << " [ " << ii << " " << jj << " 0 0 ] "<< endl;
        for ( int kk = 0 ; kk != norb ; kk++) { 
          if ( jj == kk ) {
            for ( int ll = 0 ; ll != norb ; ll++)  
               cout << gammaN[pos++] - gamma2[ii*norb+ll] << " " ;
          } else if ( jj != kk ) {
            for ( int ll = 0 ; ll != norb ; ll++) 
               cout << gammaN[pos++] <<  " " ;
          }
          cout << endl; 
        }
      }
    }

  
  } 
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::compute_sigmaN( shared_ptr<GammaInfo> gammaN_info )  {
////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::compute_sigmaN" << endl;

  string  Bra_name = gammaN_info->Bra_name();
  string  Ket_name = gammaN_info->Prev_Bra_name();

  if ( Bra_name == Ket_name ) { 

    get_wfn_data( gammaN_info->prev_Bra_info );

    shared_ptr<Determinants> Ket_det = det_old_map->at( Ket_name ); 

    int sorder  = gammaN_info->order; 
    int orb_dim = pow(Ket_det->norb(), sorder-2);
    int orb2    = Ket_det->norb()*Ket_det->norb();

    shared_ptr<Dvec> prev_sigma = dvec_sigma_map->at( gammaN_info->prev_sigma_name() );
    shared_ptr<Dvec> sigmaN     = make_shared<Dvec>( Ket_det, orb_dim*orb2  );

    for ( int  ii = 0; ii != orb_dim; ii++) {
      sigma_2a1( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), Ket_det );
      sigma_2a2( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), Ket_det );
    }

    dvec_sigma_map->emplace( gammaN_info->sigma_name, sigmaN );

    shared_ptr<Civec>  Ket = cvec_old_map->at( Ket_name );  

  } else {
  
    cout << "spin transitions sigmas not implemented yet " << endl;  
  
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::get_gamma2_from_sigma2( shared_ptr<GammaInfo> gamma2_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "B_Gamma_Computer::get_gamma2_from_sigma2_and_civec" << endl; 


  string sigma2_name =  gamma2_info->sigma_name;
  string Bra_name    =  gamma2_info->Bra_name();
  shared_ptr<Dvec>         sigma2 = dvec_sigma_map->at( sigma2_name );  
  shared_ptr<Civec>        Bra   = cvec_old_map->at( Bra_name );      
  shared_ptr<Determinants> dets   = det_old_map->at( Bra_name ); 
 
  int norb = dets->norb();

  unique_ptr<double[]> gamma2_data(new double[norb*norb]);
  for ( int ii = 0 ; ii != norb ; ii++) 
    for ( int jj = 0 ; jj != norb ; jj++) 
      *(gamma2_data.get() + (ii*norb+jj)) =  ddot_( sigma2->det()->size(),  sigma2->data(ii*norb+jj)->data(), 1, Bra->data(), 1); 

  cout << "gamma2 old method " << endl;
  for ( int ii = 0 ; ii != norb ; ii++) {
    for ( int jj = 0 ; jj != norb ; jj++) { 
      cout << *(gamma2_data.get()+ii*norb+jj) << " " ;cout.flush();
    }
    cout <<endl;
  }

////////////////////////////////////////////////////////////////////////////
  unique_ptr<double[]> gamma4_from_sigma2(new double[norb*norb*norb*norb]);
  cout << "printing gamma4 from old method sigma2 contraction ijlk" << endl;
  int pos = 0 ; 
  for ( int ii = 0 ; ii != norb ; ii++) {
    for ( int jj = 0 ; jj != norb ; jj++) { 
      cout << " [ " << ii << " " << jj << " 0 0 ] "<< endl;
        for ( int  kk = 0; kk != norb; kk++) {
          for ( int  ll = 0; ll != norb; ll++) {
            gamma4_from_sigma2[pos] = ddot_( sigma2->det()->size(), sigma2->data(ii*norb+jj)->data(),1, sigma2->data(ll*norb+kk)->data(), 1);
            cout << gamma4_from_sigma2[pos++] << " " ;
        }
        cout << endl; 
      }
    }
  }
  cout << "doing tensor contraction " << endl;
  shared_ptr<Tensor_<double>> Tens_sigma2 = Sigma_data_map->at( sigma2_name );
  convert_civec_to_tensor( Bra_name );
  shared_ptr<vector<IndexRange>> gamma_ranges = Get_Bagel_IndexRanges( gamma2_info->id_ranges ); 
  shared_ptr<Tensor_<double>> Tens_gamma2 = Tensor_Calc->contract_tensor_with_vector( Tens_sigma2, CIvec_data_map->at( Bra_name ), make_pair(gamma2_info->order, 0) );
  Gamma_data_map->emplace( gamma2_info->name, Tens_gamma2 ); 

  Print_Tensor( Tens_gamma2 , gamma2_info->name ) ;

  cout << "leaving B_Gamma_Computer::B_Gamma_Computer::get_gamma2_from_sigma2 " <<endl;
  return; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::compute_sigma2( shared_ptr<GammaInfo> gamma2_info )  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::compute_sigma2" << endl;

  string Bra_name = gamma2_info->Bra_name();
  string Ket_name = gamma2_info->Ket_name();

  get_wfn_data( gamma2_info->Ket_info );

  if ( Bra_name == Ket_name ) { 

    cout << "Ket_name = " << Ket_name << endl; 
    shared_ptr<Determinants> Ket_det = det_old_map->at( Ket_name ); 
    shared_ptr<Civec>        Ket = cvec_old_map->at( Ket_name );
    
    shared_ptr<Dvec> sigma2 = make_shared<Dvec>( Ket_det, Ket_det->norb()*Ket_det->norb() );
    
    sigma_2a1( Ket->data(), sigma2->data(0)->data(), Ket_det );
    sigma_2a2( Ket->data(), sigma2->data(0)->data(), Ket_det );
    
    dvec_sigma_map->emplace( gamma2_info->sigma_name, sigma2 );

  } else {

    cout << "spin transitions sigmas not implemented yet " << endl;  

  }
  return;
}
  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::sigma_2a1(double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets  )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                               
  //cout << "sigma_2a1" << endl;                                                                                          
  const int lb = dets->lenb();                                                                                            
  const int ij = dets->norb()*dets->norb();                                                                               
                                                                                                                          
  for (int ip = 0; ip != ij; ++ip) {                                                                                      
    double* target_base = sigma_ptr+dets->size()*ip;                                                                              
    for (auto& iter : dets->phia(ip)) {                                                                                   
      const double sign = static_cast<double>(iter.sign);                                                                 
      double* const target_array = target_base + iter.source*lb;                                                          
      blas::ax_plus_y_n(sign, cvec_ptr + iter.target*lb, lb, target_array);                                               
    }                                                                                                                     
  }                                                                                                                       
}                                                                                                                         
                                                                                                                          
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                             
void B_Gamma_Computer::B_Gamma_Computer::sigma_2a2( double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                             
//  cout << "sigma_2a2" << endl;
  const int la = dets->lena();
  const int lb = dets->lenb();
  const int ij = dets->norb()*dets->norb();

  for (int i = 0; i < la; ++i) {
    double* source_array0 = cvec_ptr+i*lb;

    for (int ip = 0; ip != ij; ++ip) {
      double* target_array0 = sigma_ptr + (ip * dets->size()) + i*lb;

      for (auto& iter : dets->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
 
                                                                                                                                                    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::get_wfn_data( shared_ptr<CIVecInfo<double>>  cvec_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string cvec_name = cvec_info->name();
  auto cvec_det_loc = det_old_map->find(cvec_name); 
  if( cvec_det_loc == det_old_map->end()){
    int II = cvec_info->state_num();
    shared_ptr<const Determinants> det_cvec_orig = cc_->data(II)->det();
    shared_ptr<Determinants> det_cvec = make_shared<Determinants>(det_cvec_orig->norb(), det_cvec_orig->nelea(), det_cvec_orig->neleb(), false, /*mute=*/true);
    det_old_map->emplace( cvec_name, det_cvec);
    cvec_old_map->emplace( cvec_name, make_shared<Civec>(*(cc_->data(II))) );
  }  

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void B_Gamma_Computer::B_Gamma_Computer::convert_civec_to_tensor( string civec_name )  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "B_Gamma_Computer::convert_civec_to_tensor" << endl;

  if ( CIvec_data_map->find(civec_name) == CIvec_data_map->end()){
 
    vector<IndexRange> civec_idxrng(1, *(range_conversion_map->at(civec_name)) );  
    
    cout <<" civec_name = " << civec_name << endl;
    cout <<" civec_idxrng[0].nblock()       = " << civec_idxrng[0].nblock()     <<  endl;
    cout <<" civec_idxrng[0].size()         = " << civec_idxrng[0].size()       <<  endl;
    cout <<" civec_idxrng[0].range().size() = " << civec_idxrng[0].range().size() <<  endl;
    cout <<" civec_idxrng[0].range().size() = " << civec_idxrng[0].range().size() <<  endl;
    
    shared_ptr<Tensor_<double>> civec_tensor = make_shared<Tensor_<double>>( civec_idxrng );
    civec_tensor->allocate();
    civec_tensor->zero();
    
    size_t idx_position = 0;
    shared_ptr<Civec> civector = cvec_old_map->at(civec_name); 
    
    for ( Index idx_block : civec_idxrng[0].range() ){
       unique_ptr<double[]> civec_block(new double[idx_block.size()]);
       std::fill_n(civec_block.get(), idx_block.size(), 0.0);
       copy_n( civector->data() + idx_position, idx_block.size(), civec_block.get());
    
       civec_tensor->put_block(civec_block, vector<Index>({ idx_block })) ;  
       idx_position += idx_block.size();  
    }
    
    CIvec_data_map->emplace( civec_name, civec_tensor); 

  }
  return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> B_Gamma_Computer::B_Gamma_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "B_Gamma_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(*ranges_str, "ranges_str" ) ; cout << endl;

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}

#endif
