#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;
#define __DEBUG_B_GAMMA_COMPUTER
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO this version is non-relativistic, must switch back to relativistic version with tensors at some point
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
B_Gamma_Computer::B_Gamma_Computer<DataType>::B_Gamma_Computer( std::shared_ptr<const Dvec > civectors ) : cc_(civectors){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::B_Gamma_Computer" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  civec_maxtile_ = 100000;// TODO set these elsewhere
  thresh_ = 1.0e-12;

  tensor_calc_ = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  dvec_sigma_map_ = make_shared<std::map< std::string, std::shared_ptr<Dvec>>>();
  det_old_map_    = make_shared<std::map< std::string, std::shared_ptr<Determinants>>>();
  cvec_old_map_   = make_shared<std::map< std::string, std::shared_ptr<Civec>>>();

  fill_and_link_determinant_map( cc_->det()->nelea() + cc_->det()->neleb(), cc_->det()->norb() );

} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::fill_and_link_determinant_map( int nelec, int norb ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::fill_and_link_determinant_map" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // if spin == 0 we are changing alphas, if spin == 1 we are changing betas
  // plusdet is the determinant with more of the kind of electron we are annihilating
  int min_spin = nelec < norb ? 0 : (nelec - norb) ; 
  int max_spin = nelec > norb ? norb : nelec; 
  
  int num_a = min_spin;
  int num_b = max_spin;

  if (min_spin != max_spin ) {
    shared_ptr<Determinants> det = make_shared<Determinants>( norb, num_a, num_b, false /*compress*/, true /*mute*/);
    for ( ;  num_a != max_spin; num_a++, num_b-- ) { 
      shared_ptr<Determinants> plus_det = make_shared<Determinants>( norb, num_a+1, num_b-1, false /*compress*/, true /*mute*/);
     
      cout << " det->nelea()  = " << det->nelea() ; cout.flush(); cout << " det->neleb()  = " << det->neleb() <<endl;
      cout << " plus_det->nelea() = " << plus_det->nelea(); cout.flush(); cout << " plus_det->neleb() = " << plus_det->neleb() << endl;

      if ( plus_det->nelea() != 0 && det->nelea() != norb )  { 

        CIStringSpace<CIStringSet<FCIString>> space({ plus_det->stringspacea(), det->stringspacea() } );
      
        space.build_linkage(1);
      
        plus_det->set_remalpha(det);
        plus_det->set_phidowna(space.phidown(plus_det->stringspacea()));
      
        det->set_addalpha(plus_det);
        det->set_phiupa(space.phiup(det->stringspacea()));

      }
      if ( det->neleb() != 0 && plus_det->neleb() != norb )  { 

        //  if the number of alpha electrons are uneven, set factor to -1; 
        const int fac = ( plus_det->nelea() & 1 ) ? -1 : 1;
        CIStringSpace<CIStringSet<FCIString>> space( { det->stringspaceb(), plus_det->stringspaceb() } );
        space.build_linkage(fac);
    
        det->set_rembeta(plus_det);
        det->set_phidownb(space.phidown(det->stringspaceb()));
      
        plus_det->set_addbeta(det);
        plus_det->set_phiupb(space.phiup(plus_det->stringspaceb()));
      }
      det_old_map_->emplace( get_det_name( 'a', det->nelea(), 'A' , det->neleb(), det->norb() ), det ); 
      det = plus_det; 
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Temporary function; ideally B_Gamma_Computer should be constructed inside System_Computer, but the dependence on the Dvecs is preventing this
// So for the time being construct B_gamma_computer with the old maps in proptool, feed this to system computer, and then alter the gamma_info_map_
// as appropriate. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map,
                                                        std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> gamma_info_map,
                                                        std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> gamma_data_map,
                                                        std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> sigma_data_map,  
                                                        std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> civec_data_map ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << " B_Gamma_Computer::B_Gamma_Computer<DataType>::Set_maps" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  range_conversion_map_ = range_conversion_map;
  gamma_info_map_ = gamma_info_map;
  gamma_data_map_ = gamma_data_map;
  sigma_data_map_ = sigma_data_map;
  civec_data_map_ = civec_data_map;
 
  new_sigma_data_map_ = make_shared<map<string, shared_ptr<Vector_Bundle<DataType>> >>();

  return;
}
/////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::get_gamma( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::get_gamma :  " <<  gamma_name << endl;
#endif //////////////////////////////////////////////////////////////////////////////

  if ( gamma_data_map_->find( gamma_name ) == gamma_data_map_->end()) {

    shared_ptr<GammaInfo_Base> gamma_info = gamma_info_map_->at(gamma_name);

    //note this has reverse iterators!
    if (gamma_info->order() > 2 ) { 
    
      get_gamma(gamma_info->prev_gamma_name());
      compute_sigmaN( gamma_info );
      convert_Dvec_sigma_to_tensor( gamma_info );
      get_gammaN_from_sigmaN ( gamma_info ) ;

    } else {

      compute_sigma2( gamma_info );
      convert_Dvec_sigma_to_tensor( gamma_info );
      get_gamma2_from_sigma2( gamma_info );

      sigma_ab( gamma_info );
      throw logic_error("got sigma_ab, now die");

    }
  }
  return;                              
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Obviously a hack for now; duplicates everything.
// Assumes that cimaxblock is _larger_ than the length of the civec. 
// With cimaxblock = 10,000 this is OK up to (8,8), so big enough for lanthanides
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::convert_Dvec_sigma_to_tensor( shared_ptr<GammaInfo_Base> gamma_info ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::Convert_to_sigma_dvec_tensor :" << gamma_info->sigma_name() <<  endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string sigma_name = gamma_info->sigma_name();
  int order = gamma_info->order();
  shared_ptr<Dvec> Dvec_sigma = dvec_sigma_map_->at( sigma_name ) ;

  shared_ptr<vector<IndexRange>> sigma_ranges = Get_Bagel_IndexRanges( gamma_info->sigma_id_ranges() ); 

  shared_ptr<Tensor_<DataType>> Tens_sigma = make_shared<Tensor_<DataType>>( *sigma_ranges ); 
  Tens_sigma->allocate();
  Tens_sigma->zero();
  sigma_data_map_->emplace( sigma_name, Tens_sigma ); 

  shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( *sigma_ranges ) ;

  vector<int> range_lengths = get_range_lengths( *sigma_ranges ) ;
  vector<int> block_pos(range_lengths.size(),0);  
  vector<int> mins(range_lengths.size(),0);  

  vector<int> orb_id_strides(order);
  for ( int ii = 0 ; ii != order ; ii++ ) 
    orb_id_strides[ii] = (int)pow( Dvec_sigma->det()->norb(), ii );


  auto boffset = [&block_offsets, &block_pos]( int pos ){ return block_offsets->at(pos).at(block_pos.at(pos)); }; 
  do {
      
    vector<Index> sigma_id_blocks = get_rng_blocks( block_pos, *sigma_ranges );
    int block_size = sigma_id_blocks[0].size();
    int orb_pos = 0;
    //TODO check: starts at one because ci_index is the first index, although do not understand  how this was working before
    //            CHECK THIS; it's a change made upon returning to this section of the code after work elsewhere.
    for ( int  ii = 1 ; ii != order+1; ii++ ){
      orb_pos += boffset(ii)*orb_id_strides[ii-1]; 
      block_size *= sigma_id_blocks[ii].size(); 
    }    
    unique_ptr<DataType[]> sigma_block( new DataType[block_size] );
    fill_n( sigma_block.get(), block_size, DataType(0.0) );
    DataType* sigma_block_ptr = sigma_block.get();

    int ci_block_size = sigma_id_blocks[0].size();
    int sblock_pos = 0; 
  
    //TODO almost certain boffset(order) is wrong, should be block_offset(0), needs to be offset of ci index 
    int num_orb_idxs_in_this_block = block_size/ci_block_size;
    for (int ii = 0; ii != num_orb_idxs_in_this_block; ii++ ) { 
      std::copy( Dvec_sigma->data(orb_pos)->data()+boffset(0), Dvec_sigma->data(orb_pos)->data()+boffset(0)+ci_block_size, sigma_block_ptr );   
      sigma_block_ptr+=ci_block_size;
      orb_pos++;
    }

    Tens_sigma->put_block( sigma_block, sigma_id_blocks ) ; 

  } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::get_gammaN_from_sigmaN( shared_ptr<GammaInfo_Base> gammaN_info ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::get_gammaN_from_sigmaN : " << gammaN_info->name() <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  string gammaN_name = gammaN_info->name();
  shared_ptr<Dvec> sigmaN = dvec_sigma_map_->at(gammaN_info->sigma_name());
  string Bra_name = gammaN_info->Bra_name();
  shared_ptr<Civec> Bra = cvec_old_map_->at( Bra_name );      
  shared_ptr<Determinants> Bra_det = det_old_map_->at( Bra_name ) ; 

  int orb_dim = pow(Bra_det->norb(), (gammaN_info->order() - 2) );
  int norb    = Bra_det->norb();
  int orb2    = norb*norb;
  int order = gammaN_info->order();

  shared_ptr<vector<IndexRange>> gammaN_ranges = Get_Bagel_IndexRanges( gammaN_info->id_ranges() ); 
  shared_ptr<Tensor_<DataType>> Tens_gammaN = make_shared<Tensor_<DataType>>( *gammaN_ranges ); 
  Tens_gammaN->allocate();
  Tens_gammaN->zero();
  gamma_data_map_->emplace( gammaN_name, Tens_gammaN ); 

  shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( *gammaN_ranges ) ;

  vector<int> range_lengths  =  get_range_lengths( *gammaN_ranges );
  vector<int> block_pos(order,0);  
  vector<int> mins(order,0);  

  vector<int> orb_id_strides(order);
  for ( int ii = 0 ; ii != order ; ii++ ) 
    orb_id_strides[ii] = (int)pow(norb, ii);

  auto boffset = [&block_offsets, &block_pos]( int pos ){ return block_offsets->at(pos).at(block_pos.at(pos)); }; 
  do {
      
    vector<Index> gammaN_id_blocks = get_rng_blocks( block_pos, *gammaN_ranges );
    int block_size = 1 ;
    int orb_pos = 0;
    for ( int  ii = 0 ; ii != order; ii++ ){
      orb_pos += boffset(ii)*orb_id_strides[ii]; 
      block_size *= gammaN_id_blocks[ii].size(); 
    }    

    unique_ptr<DataType[]> gammaN_block( new DataType[block_size] );
    fill_n( gammaN_block.get(), block_size, DataType(0.0) );
    DataType* gammaN_block_ptr = gammaN_block.get();

    for (int ii = orb_pos; ii != orb_pos+block_size; ii++ )  
      *gammaN_block_ptr++ = ddot_( Bra_det->size(),sigmaN->data(ii)->data(), 1, Bra->data(), 1); 

    Tens_gammaN->put_block( gammaN_block, gammaN_id_blocks ) ; 

  } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This assumes all orbitals have the same size of range, OK for alpha beta spins, not OK 
// if different active spaces.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::compute_sigmaN( shared_ptr<GammaInfo_Base> gammaN_info )  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::compute_sigmaN : " << gammaN_info->sigma_name() ;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gammaN_info->Bra_name();     
  string ket_name = gammaN_info->Prev_Bra_name();
  get_wfn_data( gammaN_info->prev_Bra_info() );
  shared_ptr<Determinants> ket_det = det_old_map_->at( ket_name );  

  int sorder  = gammaN_info->order(); 
  int orb_dim = pow(ket_det->norb(), sorder-2);
  int orb2    = ket_det->norb()*ket_det->norb();

  if ( dvec_sigma_map_->find(gammaN_info->prev_sigma_name()) == dvec_sigma_map_->end() ) {
    if ( sorder > 4 ) {
      compute_sigmaN( gamma_info_map_->at(gammaN_info->prev_gamma_name()) ); 
    } else { 
      compute_sigma2( gamma_info_map_->at(gammaN_info->prev_gamma_name()) ); 
    }
  } 
  shared_ptr<Dvec> prev_sigma = dvec_sigma_map_->at( gammaN_info->prev_sigma_name() ); 
  shared_ptr<Dvec> sigmaN     = make_shared<Dvec>( ket_det, orb_dim*orb2  );

  if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha() ){ 

    for ( int  ii = 0; ii != orb_dim; ii++) {
      sigma_2a1( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), ket_det );
      sigma_2a2( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), ket_det );
    }

    dvec_sigma_map_->emplace( gammaN_info->sigma_name(), sigmaN );

  } else if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()+1 ){ 

    shared_ptr<Determinants> bra_det = det_old_map_->at( bra_name );  
    for ( int  ii = 0; ii != orb_dim; ii++) 
      sigma_2ab( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), bra_det, ket_det );
  
  } else if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()-1 ){ 
    
    throw logic_error( " sigma for (+a,-b) not implemented yet " ) ; 

  } else { 
    cout << "The density matrix transition required by the following is not implemented " << endl;
    cout << "gammaN_info->Bra_nalpha()       = "; cout.flush(); cout << gammaN_info->Bra_nalpha() << endl;
    cout << "gammaN_info->prev_Bra_nalpha()  = "; cout.flush(); cout << gammaN_info->prev_Bra_nalpha()<< endl; 
    cout << "gammaN_info->Bra_nbeta()       = "; cout.flush(); cout << gammaN_info->Bra_nbeta() << endl;
    cout << "gammaN_info->prev_Bra_nbeta()  = "; cout.flush(); cout << gammaN_info->prev_Bra_nbeta()<< endl; 
    throw logic_error(" aborting" );
  }
  
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::get_gamma2_from_sigma2( shared_ptr<GammaInfo_Base> gamma2_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::get_gamma2_from_sigma2" << endl; 
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string sigma2_name =  gamma2_info->sigma_name();
  string bra_name    =  gamma2_info->Bra_name();
 
  convert_civec_to_tensor( bra_name );

  shared_ptr<Tensor_<DataType>> gamma2_tens = tensor_calc_->contract_tensor_with_vector( sigma_data_map_->at( sigma2_name ), civec_data_map_->at( bra_name ), 0 );
  gamma_data_map_->emplace( gamma2_info->name(), gamma2_tens ); 

  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::compute_sigma2_vb( shared_ptr<GammaInfo_Base> gamma2_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::compute_sigma2_vb : " << gamma2_info->name() << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  string ket_name = gamma2_info->Ket_name();
  get_wfn_data( gamma2_info->Ket_info() );

  shared_ptr<Determinants> ket_det = det_old_map_->at( ket_name ); 
  shared_ptr<Civec>        ket = cvec_old_map_->at( ket_name );

  convert_civec_to_tensor( ket_name );

  if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha() ) { 

   sigma_aa_vb( gamma2_info, true  );
   sigma_bb_vb( gamma2_info, false );
    
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()+1 ) { 
   sigma_ab_vb( gamma2_info, true );

  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()-1 ) { 
   sigma_ba_vb( gamma2_info, true );

  } else {
    throw logic_error( "this sigma: " + gamma2_info->sigma_name() + " is not implemented" ); 
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::compute_sigma2( shared_ptr<GammaInfo_Base> gamma2_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::compute_sigma2 : " << gamma2_info->name() << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  string ket_name = gamma2_info->Ket_name();

  get_wfn_data( gamma2_info->Ket_info() );

  shared_ptr<Determinants> ket_det = det_old_map_->at( ket_name ); 
  shared_ptr<Civec>        ket = cvec_old_map_->at( ket_name );
  
  shared_ptr<Dvec> sigma2;
 
  if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha() ) { 

    sigma2 = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );
    sigma_2a1( ket->data(), sigma2->data(0)->data(), ket_det );
    sigma_2a2( ket->data(), sigma2->data(0)->data(), ket_det );

    sigma2_test( gamma2_info );
   
    
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()+1 ) { 

     shared_ptr<Determinants> bra_det = det_old_map_->at( gamma2_info->Bra_name() );

     sigma2 = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );

     sigma_2ab( ket->data(), sigma2->data(0)->data(), bra_det, ket_det );
  
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()-1 ) { 

    throw logic_error( "(-a,+b) sigmas not implemented yet. Aborting !!  ");  

  } else {

    cout << "this sigma : "; cout.flush(); cout << gamma2_info->sigma_name() << " is not implemented" << endl;
    throw logic_error( "Aborting !!  ");  

  }

  dvec_sigma_map_->emplace( gamma2_info->sigma_name(), sigma2 );

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2a1(DataType* cvec_ptr, DataType* sigma_ptr, shared_ptr<Determinants> dets  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2a1" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const int lb = dets->lenb();
  const int ij = dets->norb()*dets->norb();
  for (int ip = 0; ip != ij; ++ip) {
    DataType* target_base = sigma_ptr+dets->size()*ip;
    for (auto& iter : dets->phia(ip)) {
      const DataType sign = static_cast<DataType>(iter.sign);
      DataType* const target_array = target_base + iter.source*lb;
      blas::ax_plus_y_n(sign, cvec_ptr + iter.target*lb, lb, target_array);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2a2( DataType* cvec_ptr, DataType* sigma_ptr, shared_ptr<Determinants> dets) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2a2" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  const int la = dets->lena();
  const int lb = dets->lenb();
  const int ij = dets->norb()*dets->norb();
  for (int i = 0; i < la; ++i) {
    DataType* source_array0 = cvec_ptr+i*lb;
    for (int ip = 0; ip != ij; ++ip) {
      DataType* target_array0 = sigma_ptr + (ip * dets->size()) + i*lb;
      for (auto& iter : dets->phib(ip)) {
        const DataType sign = static_cast<DataType>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab( shared_ptr<GammaInfo_Base> gamma_info ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  shared_ptr<Determinants>  bra_det = det_old_map_->at(bra_name);

  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants>  ket_det = det_old_map_->at(ket_name);
  bra_det = make_shared<Determinants>(ket_det->norb(), ket_det->nelea()+1, ket_det->neleb()-1, false, false);

  assert(ket_det->nelea() == bra_det->nelea()-1 && ket_det->neleb()== bra_det->neleb()+1 );

  auto ket = cvec_old_map_->find(ket_name); 

  const int norb = ket_det->norb();
  const int ket_lena = ket_det->lena();
  const int ket_lenb = ket_det->lenb();
  const int bra_lena = bra_det->lena();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  shared_ptr<Dvec> sigma_ab = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  DataType* ket_ptr_buff = cvec_old_map_->at(ket_name)->data();
  for (int ii = 0; ii < norb; ++ii) {
    for (int jj = 0; jj < norb; ++jj) {
      unique_ptr<DataType[]> bra_data( new DataType[bra_length] );
      fill_n( bra_data.get() , bra_length, (DataType)(0.0) );
      DataType* bra_ptr = bra_data.get();
      DataType* ket_ptr = ket_ptr_buff;
      for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
        for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it, ++ket_ptr ) {
          if ( !(*abit_it)[ii] && (*bbit_it)[jj] ) {
            bitset<nbit__> abit_bra = *abit_it;
            bitset<nbit__> bbit_bra = *bbit_it;
            abit_bra.set(ii);
            bbit_bra.reset(jj);
            const double abphase =  ket_det->sign<0>( (*abit_it), ii) * ket_det->sign<1>( (*bbit_it), jj);
            *(bra_ptr + bra_det->lexical<0>(abit_bra)*bra_lenb + bra_det->lexical<1>(bbit_bra)) += abphase * (*ket_ptr);
          }
        }
      }
      move( bra_ptr, bra_ptr+bra_length, sigma_ab->data(ii*norb + jj)->data() );
    }
  }

//  dvec_sigma_map_->emplace( gamma_info->name(), sigma_ab );
  
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_aa_test( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_aa_test" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = ket_det;

  auto ket = cvec_old_map_->find(ket_name); 

  const int norb = ket_det->norb();
  const int ket_lena = ket_det->lena();
  const int ket_lenb = ket_det->lenb();
  const int bra_lena = bra_det->lena();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  shared_ptr<Dvec> sigma_aa;
  if ( new_sigma) {
    sigma_aa = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  } else { 
    sigma_aa = dvec_sigma_map_->at( gamma_info->sigma_name()) ;
  }

  DataType* ket_ptr_buff = cvec_old_map_->at(ket_name)->data();

  for (int ii = 0; ii != norb; ii++) {
    for (int jj = 0; jj != norb; jj++) {
      for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
        bool possible_exc = ( jj == ii ) ? (*abit_it)[jj] : !(*abit_it)[ii] && (*abit_it)[jj]; 
        if ( possible_exc ) {
            
          DataType* aiaj_ket_section_start = ket_ptr_buff + ket_det->lexical<0>( *abit_it ) * ket_lenb;

          bitset<nbit__> abit_bra = *abit_it;
          DataType op_phase =  ket_det->sign<0>( abit_bra, jj);
          abit_bra.reset(jj);
          op_phase *= ket_det->sign<0>( abit_bra, ii);
          abit_bra.set(ii);

          DataType* bra_section_start = sigma_aa->data(ii + norb*jj)->data() +  bra_det->lexical<0>( abit_bra ) * bra_lenb;
          blas::ax_plus_y_n( op_phase, aiaj_ket_section_start, bra_lenb, bra_section_start);
        }
      }
    }
  }

  dvec_sigma_map_->emplace( "sigma_aa_test", sigma_aa );
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_bb_test( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_bb_test" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = ket_det;

  auto ket = cvec_old_map_->find(ket_name); 

  const int norb = ket_det->norb();
  const int ket_lena = ket_det->lena();
  const int ket_lenb = ket_det->lenb();
  const int bra_lena = bra_det->lena();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  shared_ptr<Dvec> sigma_bb = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  if ( new_sigma) {
    sigma_bb = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  } else { 
    sigma_bb = dvec_sigma_map_->at( "sigma_bb_test" ) ;
  }

//  DataType trans_factor = (ket_det->nelea()*ket_det->neleb() & 1)  ? (DataType)(-1.0) : (DataType)(1.0);
  unique_ptr<DataType[]> ket_transposed_data( new DataType[ket_lena*ket_lenb] );
  fill_n( ket_transposed_data.get(), ket_lena*ket_lenb, (DataType)(0.0) );
  DataType* ket_ptr_buff = ket_transposed_data.get();
  blas::transpose( cvec_old_map_->at(ket_name)->data() , ket_lenb, ket_lena, ket_ptr_buff);

  for (int ii = 0; ii != norb; ii++) {
    for (int jj = 0; jj != norb; jj++) {
      unique_ptr<DataType[]> bibj_ket_data( new DataType[bra_length] );
      fill_n(bibj_ket_data.get(), bra_length, (DataType)(0.0));
      DataType* bb_ket_ptr = bibj_ket_data.get();
      
      for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it ) { // loop through a strings
        bool possible_exc = ( jj == ii ) ? (*bbit_it)[jj] : !(*bbit_it)[ii] && (*bbit_it)[jj]; 
        if ( possible_exc ) {
            
          DataType* bibj_ket_section_start = ket_ptr_buff + ket_det->lexical<0>( *bbit_it ) * ket_lena;

          bitset<nbit__> bbit_bra = *bbit_it;                   
          DataType op_phase =  ket_det->sign<0>( bbit_bra, jj); // using sign<0> instead of sign<1> to avoid scaling of ket_data
          bbit_bra.reset(jj);
          op_phase *= ket_det->sign<0>( bbit_bra, ii);
          bbit_bra.set(ii);

          DataType* bra_section_start = bb_ket_ptr + bra_det->lexical<0>( bbit_bra ) * bra_lena;  // positions of this section in the bra
          blas::ax_plus_y_n( op_phase, bibj_ket_section_start, bra_lena, bra_section_start);
        }
      }
      blas::transpose(bibj_ket_data.get(), bra_lena, bra_lenb, sigma_bb->data(ii+norb*jj)->data() );
    }
  }

  dvec_sigma_map_->emplace( "sigma_bb_test", sigma_bb );
  
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab_test( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab_test" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = det_old_map_->at(bra_name);

  bra_det = make_shared<Determinants>( ket_det->norb(), ket_det->nelea()+1, ket_det->neleb()-1, false /*compress*/, true /*mute*/);
  assert ( bra_det->nelea() == ket_det->nelea()+1  && bra_det->neleb() == ket_det->neleb()-1 );

  auto ket = cvec_old_map_->at(ket_name); 
  DataType* ket_ptr_orig = ket->data();

  const int norb = ket_det->norb();
  const int ket_lenb = ket_det->lenb();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  shared_ptr<Dvec> sigma_ab;
  if ( new_sigma) {
    sigma_ab = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  } else { 
    sigma_ab = dvec_sigma_map_->at( gamma_info->sigma_name() ) ;
  }

  // Slow, but simple to check and parallelize
  for (int ii = 0; ii != norb; ii++) {
    for (int jj = 0; jj != norb; jj++) {
      DataType* sigma_ab_ij_ptr = sigma_ab->data(ii +norb*jj)->data();
      DataType* ket_ptr = ket_ptr_orig; 
      for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
        for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it, ++ket_ptr ) {
          if ( !(*abit_it)[ii] && (*bbit_it)[jj] ) {

            bitset<nbit__> bbit_bra = *bbit_it;
            DataType op_phase =  ket_det->sign<1>( bbit_bra, jj);
            bbit_bra.reset(jj);

            bitset<nbit__> abit_bra = *abit_it;
            op_phase *= ket_det->sign<0>( abit_bra, ii);
            abit_bra.set(ii);
            { 
              DataType* bra_ptr = sigma_ab_ij_ptr + bra_det->lexical<0>( abit_bra ) * bra_lenb + bra_det->lexical<1>( bbit_bra ); 
              *bra_ptr += op_phase* (*ket_ptr);
            }
          }
        }
      }
    }
  }

  dvec_sigma_map_->emplace( "sigma_ab_test", sigma_ab );
  
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ba_test( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ba_test" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = det_old_map_->at(bra_name);

  bra_det = make_shared<Determinants>( ket_det->norb(), ket_det->nelea()-1, ket_det->neleb()+1, false /*compress*/, true /*mute*/);
  assert ( bra_det->nelea() == ket_det->nelea()-1  && bra_det->neleb() == ket_det->neleb()+1 );

  auto ket = cvec_old_map_->at(ket_name); 
  DataType* ket_ptr_orig = ket->data();

  const int norb = ket_det->norb();
  const int ket_lenb = ket_det->lenb();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  shared_ptr<Dvec> sigma_ba;
  if ( new_sigma) {
    sigma_ba = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );
  } else { 
//    sigma_ba = dvec_sigma_map_->at( gamma_info->sigma_name()) ;
  }

  // Slow, but simple to check and parallelize
  for (int ii = 0; ii != norb; ii++) {
    for (int jj = 0; jj != norb; jj++) {
      DataType* sigma_ba_ij_ptr = sigma_ba->data(ii +norb*jj)->data();
      DataType* ket_ptr = ket_ptr_orig; 
      for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
        for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it, ++ket_ptr ) {
          if (  !(*bbit_it)[ii] && (*abit_it)[jj]) {

            bitset<nbit__> abit_bra = *abit_it;
            DataType op_phase = ket_det->sign<0>( abit_bra, jj);
            abit_bra.reset(jj);

            bitset<nbit__> bbit_bra = *bbit_it;
            op_phase *= ket_det->sign<1>( bbit_bra, ii);
            bbit_bra.set(ii);

            { 
              DataType* bra_ptr = sigma_ba_ij_ptr + bra_det->lexical<0>( abit_bra ) * bra_lenb + bra_det->lexical<1>( bbit_bra ); 
              *bra_ptr += op_phase* (*ket_ptr);
            }
          }
        }
      }
    }
  }

  dvec_sigma_map_->emplace( gamma_info->sigma_name(), sigma_ba );
  dvec_sigma_map_->emplace( "sigma_ba_test", sigma_ba );
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2ab( DataType* cvec_ptr, DataType* sigma_ptr, 
                                                              shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2ab" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //get pointer in target tensor
  size_t bra_civec_length = bra_det->lena()* bra_det->lenb();
  size_t bra_lenb = bra_det->lenb();
  size_t ket_lenb = ket_det->lenb();
  size_t norb = bra_det->norb();
  for (unsigned int ii =0 ; ii!=norb ; ii++) { 
    for (unsigned int jj =0 ; jj!=norb ; jj++) {
      DataType* sigma_ij_0_ptr = sigma_ptr + (ii*norb+jj)*bra_civec_length;

      //get position in bra vec by doing the opposite of the operation you are going to do on the ket
      for ( auto& a1 : bra_det->phidowna(ii)){
        for ( auto& b1 : bra_det->phiupb(jj)  ){
          DataType* sigma_aj_bi_ptr = sigma_ij_0_ptr + (a1.source*bra_lenb + b1.source);

          //get position in ket by standard approach
          for ( auto& a2 : ket_det->phiupa(ii)){     
            for ( auto& b2 : ket_det->phidownb(jj)) { 
              *sigma_aj_bi_ptr += *(cvec_ptr + a2.source*ket_lenb + b2.source);
            }
          }
        }
      }
    }
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::get_wfn_data( shared_ptr<CIVecInfo_Base>  cvec_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << " B_Gamma_Computer::get_wfn_data " << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  string cvec_name = cvec_info->name();
  string det_name = cvec_info->det_name();
  auto cvec_det_loc = det_old_map_->find(cvec_name); 
  if( cvec_det_loc == det_old_map_->end()){
    det_old_map_->emplace( cvec_name, det_old_map_->at( det_name ) );
    cvec_old_map_->emplace( cvec_name, make_shared<Civec>(*(cc_->data( cvec_info->state_num() ))) );
  }  
  cout << "got " << cvec_name << " in cvec_old_map_ " <<endl; 
  cout << "got " << det_name << " in det_old_map " <<endl; 
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::convert_civec_to_tensor( string civec_name )  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::convert_civec_to_tensor" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( civec_data_map_->find(civec_name) == civec_data_map_->end()){
 
    vector<IndexRange> civec_idxrng(1, *(range_conversion_map_->at(civec_name)) );

    shared_ptr<Tensor_<DataType>> civec_tensor = make_shared<Tensor_<DataType>>( civec_idxrng );
    civec_tensor->allocate();
    civec_tensor->zero();
    size_t idx_position = 0;
    shared_ptr<Civec> civector = cvec_old_map_->at(civec_name);
    for ( Index idx_block : civec_idxrng[0].range() ){
       unique_ptr<DataType[]> civec_block(new DataType[idx_block.size()]);
       std::fill_n(civec_block.get(), idx_block.size(), DataType(0.0) );
       copy_n( civector->data() + idx_position, idx_block.size(), civec_block.get());
       civec_tensor->put_block(civec_block, vector<Index>({ idx_block })) ;
       idx_position += idx_block.size();
    }
    civec_data_map_->emplace( civec_name, civec_tensor);
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::compute_sigmaN_vb( shared_ptr<GammaInfo_Base> gammaN_info )  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::compute_sigmaN_vb : " << gammaN_info->sigma_name() ;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////

  int sorder  = gammaN_info->order(); 
  // check to see if previous sigma has been calculated, if not calculate it (recursive call here) 
  if ( new_sigma_data_map_->find(gammaN_info->prev_sigma_name()) == new_sigma_data_map_->end() ) {
    if ( sorder > 4 ) {
      compute_sigmaN( gamma_info_map_->at(gammaN_info->prev_gamma_name()) ); 
    } else { 
      compute_sigma2( gamma_info_map_->at(gammaN_info->prev_gamma_name()) ); 
    }
  }

  shared_ptr<Determinants> ket_det = det_old_map_->at( gammaN_info->Bra_name() );  
  shared_ptr<Determinants> bra_det = det_old_map_->at( gammaN_info->Prev_Bra_name() );  

  size_t bra_length = bra_det->lena()*bra_det->lenb();
  int norb = ket_det->norb();
 
  vector<int> orb_ranges_sigma_n(sorder, norb );
  shared_ptr<Vector_Bundle<DataType>> sigma_n = make_shared<Vector_Bundle<DataType>>( orb_ranges_sigma_n, bra_length, civec_maxtile_, true, false, false );

  // TODO  should define norb so can be variable from gamma_info..
  vector<int> maxs_prev_sigma(sorder-2, norb);
  vector<int> mins_prev_sigma(sorder-2, 0);
  vector<int> orb_ids_prev_sigma = mins_prev_sigma;
   

  if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha() ){ 
    vector<bool> sigma_overwrite_pattern(sorder,false);
    sigma_overwrite_pattern[0] = true;//lazy...
    sigma_overwrite_pattern[1] = true; 

    do {  
    
      vector<int> maxs_sigma2(2, norb);
      vector<int> mins_sigma2(2, 0);
      vector<int> orb_ids_sigma2 = mins_prev_sigma;
     
      do {  
         auto tmp_sigma = make_shared<Vector_Bundle<DataType>>( orb_ids_sigma2, bra_length, civec_maxtile_, true, false, false );
          
         compute_eiej_on_ket( tmp_sigma, sigma_n->vector_map(orb_ids_prev_sigma), bra_det, ket_det, "AA" ); 
         compute_eiej_on_ket( tmp_sigma, sigma_n->vector_map(orb_ids_prev_sigma), bra_det, ket_det, "BB" ); 

         sigma_n-> merge_fixed_ids( tmp_sigma, orb_ids_prev_sigma, sigma_overwrite_pattern, true );

      } while(fvec_cycle_skipper( orb_ids_sigma2, maxs_sigma2, mins_sigma2) );

    } while(fvec_cycle_skipper( orb_ids_prev_sigma, maxs_prev_sigma, mins_prev_sigma) );

  } else if ( ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()+1 ) &&
              ( gammaN_info->Bra_nbeta() ==  gammaN_info->prev_Bra_nbeta()-1 ) ){ 

    cout << "hello" << endl;
    //sigma_vb ab

  } else if ( ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()-1 ) &&
              ( gammaN_info->Bra_nbeta()  ==  gammaN_info->prev_Bra_nbeta()+1 ) ){ 
    
    cout << "hello" << endl;
    //sigma_vb ba

  } else { 
    throw logic_error( gammaN_info->prev_gamma_name() + " is not yet implemented! Aborting!!" );
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::compute_eiej_on_ket( shared_ptr<Vector_Bundle<DataType>> eiej_on_ket,
                                                                   shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                                   shared_ptr<Determinants> bra_det,
                                                                   shared_ptr<Determinants> ket_det,
                                                                   string transition_name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::compute_eiej_on_ket ";
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( transition_name == "AA") {
    assert( (bra_det->nelea() == ket_det->nelea()) && ( bra_det->neleb() == ket_det->neleb()));
    sigma2_aa_vb( eiej_on_ket, ket_tensor, bra_det, ket_det );
  
  } else if ( transition_name == "BB") {
    assert( (bra_det->nelea() == ket_det->nelea()) && ( bra_det->neleb() == ket_det->neleb()));
    sigma2_bb_vb( eiej_on_ket, ket_tensor, bra_det, ket_det );
    
  } else if ( transition_name == "AB") {
    assert( (bra_det->nelea()+1 == ket_det->nelea()) && ( bra_det->neleb()-1 == ket_det->neleb()));
    sigma2_ab_vb( eiej_on_ket, ket_tensor, bra_det, ket_det );

  } else if ( transition_name == "BA") {
    assert( (bra_det->nelea()-1 == ket_det->nelea()) && ( bra_det->neleb()+1 == ket_det->neleb()));
    sigma2_ba_vb( eiej_on_ket, ket_tensor, bra_det, ket_det );

  } else {
    throw logic_error( "B_Gamma_Computer::compute_eiej_on_ket : Aborted as this sigma is not implemented; not aa , bb, ab or ba" ); 
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_aa_vb( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_aa_test_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = det_old_map_->at(bra_name);

  shared_ptr<Vector_Bundle<DataType>> sigma_aa;
  if ( new_sigma ) {
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_aa = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );
  } else {
    sigma_aa = new_sigma_data_map_->at( gamma_info->sigma_name() );
  }
  convert_civec_to_tensor( ket_name );
  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name);
  
  sigma2_aa_vb( sigma_aa, ket_tensor, bra_det, ket_det );
  new_sigma_data_map_->emplace( "sigma_aa_vb_test", sigma_aa );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_aa_vb( shared_ptr<Vector_Bundle<DataType>> sigma_aa,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_aa_vb" << endl;
//TODO : ket_det and bra_det should be the same, only keep for now to facilitate tests on gamma task list generation
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  const int norb = ket_det->norb();
  const int ket_lenb = ket_det->lenb();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();
  
  SMITH::IndexRange ket_idx_range =  ket_tensor->indexrange().front();

  assert ( ket_idx_range.range().size() == 1 ) ; // TODO : will not work if ci block is split
  for ( SMITH::Index block : ket_idx_range.range() ) {
    vector<SMITH::Index> ket_block_id = { block };
    unique_ptr<DataType[]> ket_data = ket_tensor->get_block( ket_block_id ) ;
  
    for (int ii = 0; ii != norb; ii++) {
      for (int jj = 0; jj != norb; jj++) {
        
        vector<int> orb_ids = { ii, jj }; 
        std::shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_aa->get_new_vector( true );

        assert ( sigma_aa->index_range_vec_[0].range().size() == 1 ) ; //TODO : will not work if ci block is split

        bool ij_vec_sparse = true;
        for ( SMITH::Index& ij_block_id : sigma_aa->index_range_vec_[0] ) {
          vector<SMITH::Index> sigma_ij_block_id = { ij_block_id };
          unique_ptr<DataType[]> sigma_ij_block_data = sigma_ij_vec->get_block( sigma_ij_block_id );

          for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
            bool possible_exc = ( jj == ii ) ? (*abit_it)[jj] : !(*abit_it)[ii] && (*abit_it)[jj]; 
            if ( possible_exc ) {
          
              bitset<nbit__> abit_bra = *abit_it;
              DataType op_phase =  ket_det->sign<0>( abit_bra, jj);
              abit_bra.reset(jj);
              op_phase *= ket_det->sign<0>( abit_bra, ii);
              abit_bra.set(ii);
          
              DataType* aiaj_ket_section_start = ket_data.get() + ket_det->lexical<0>( *abit_it ) * ket_lenb;
              DataType* bra_section_start = sigma_ij_block_data.get() +  bra_det->lexical<0>( abit_bra ) * bra_lenb;

              blas::ax_plus_y_n( op_phase, aiaj_ket_section_start, bra_lenb, bra_section_start);
            }
          }
          if ( ( abs(*(max_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ ) ||
               ( abs(*(min_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ )    ) { 
            ij_vec_sparse = false;
            sigma_ij_vec->put_block( sigma_ij_block_data, sigma_ij_block_id ); 
          }
        }
        sigma_aa->set_sparsity( orb_ids, ij_vec_sparse );
        if ( !ij_vec_sparse )
          sigma_aa->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
      }
    }
    ket_tensor->put_block( ket_data );
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_bb_vb( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_bb_test_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = ket_det;

  shared_ptr<Vector_Bundle<DataType>> sigma_bb;
  if ( new_sigma ) {
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_bb = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );
  } else {
    sigma_bb = new_sigma_data_map_->at( gamma_info->sigma_name() );
  }

  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name);
  
  sigma2_bb_vb( sigma_bb, ket_tensor, bra_det, ket_det ) ;
  new_sigma_data_map_->emplace( "sigma_bb_vb_test", sigma_bb );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_bb_vb( shared_ptr<Vector_Bundle<DataType>> sigma_bb,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_bb_vb" << endl;
//TODO : ket_det and bra_det should be the same, only keep for now to facilitate tests on gamma task list generation
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  const int norb = ket_det->norb();
  const int ket_lena = ket_det->lena();
  const int ket_lenb = ket_det->lenb();
  const int bra_lena = bra_det->lena();
  const int bra_lenb = bra_det->lenb();

  vector<SMITH::IndexRange> ket_idx_range =  ket_tensor->indexrange();
  assert ( ket_idx_range.front().range().size() == 1 ) ; // TODO : will not work if ci block is split
  for ( SMITH::Index block : ket_idx_range.front().range() ) {
    vector<SMITH::Index> ket_block_id = { block };
    unique_ptr<DataType[]> ket_data_transposed( new DataType[block.size()]);
    fill_n( ket_data_transposed.get(), block.size(), (DataType)(0.0) );

    {
      unique_ptr<DataType[]> ket_data = ket_tensor->get_block( ket_block_id ) ;
      blas::transpose( ket_data.get(), ket_lenb, ket_lena, ket_data_transposed.get() );
      ket_tensor->put_block( ket_data, ket_block_id ) ;
    }

    for (int ii = 0; ii != norb; ii++) {
      for (int jj = 0; jj != norb; jj++) {
        
        vector<int> orb_ids = { ii, jj }; 
        shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_bb->get_new_vector(true);  

        assert ( sigma_bb->index_range_vec_[0].range().size() == 1 ) ; //TODO : will not work if ci block is split

        bool ij_vec_sparse = true;
        for ( SMITH::Index& ij_block_id : sigma_bb->index_range_vec_[0] ) {

          unique_ptr<DataType[]> ij_block_t( new DataType[ij_block_id.size()] ); 
          fill_n( ij_block_t.get(), ij_block_id.size(), (DataType)(0.0) ); 

          for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it ) {
            bool possible_exc = ( jj == ii ) ? (*bbit_it)[jj] : !(*bbit_it)[ii] && (*bbit_it)[jj]; 
            if ( possible_exc ) {
                
          
              bitset<nbit__> bbit_bra = *bbit_it;
              DataType op_phase =  ket_det->sign<0>( bbit_bra, jj);
              bbit_bra.reset(jj);
              op_phase *= ket_det->sign<0>( bbit_bra, ii);
              bbit_bra.set(ii);
          
              DataType* bibj_ket_section_start = ket_data_transposed.get() + ket_det->lexical<0>( *bbit_it ) * ket_lena;
              DataType* bra_section_start = ij_block_t.get() +  bra_det->lexical<0>( bbit_bra ) * bra_lena;
              blas::ax_plus_y_n( op_phase, bibj_ket_section_start, bra_lena, bra_section_start);

            }
          }

          if ( (abs(*(max_element(ij_block_t.get() , ij_block_t.get() + ij_block_id.size() ) )) > thresh_ ) || 
               (abs(*(min_element(ij_block_t.get() , ij_block_t.get() + ij_block_id.size() ))) > thresh_ )     ) { 

            ij_vec_sparse = false;
            unique_ptr<DataType[]> ij_block( new DataType[ij_block_id.size()] );
            fill_n( ij_block.get(), ij_block_id.size(), (DataType)(0.0) );
            blas::transpose( ij_block_t.get(), bra_lena, bra_lenb, ij_block.get() ); //TODO : will not work if ci block is split
            vector<SMITH::Index> ij_block_id_vec ={ ij_block_id };
            sigma_ij_vec->put_block( ij_block, ij_block_id_vec ); 
          }

        }
        sigma_bb->set_sparsity( orb_ids, ij_vec_sparse );
        if ( !ij_vec_sparse )
          sigma_bb->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
      }
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab_vb( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ab_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = det_old_map_->at(bra_name);

  {//TEST  
  bra_det = make_shared<Determinants>( ket_det->norb(), ket_det->nelea()+1, ket_det->neleb()-1, false /*compress*/, true /*mute*/);
  } // END TEST

  convert_civec_to_tensor( ket_name ); 
  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name); 

  shared_ptr<Vector_Bundle<DataType>> sigma_ab;
  if ( new_sigma) {
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_ab = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );
    new_sigma_data_map_->emplace( "sigma_ab_vb", sigma_ab );
  } else { 
    sigma_ab = new_sigma_data_map_->at( gamma_info->sigma_name() ) ;
  }
  
  sigma2_ab_vb( sigma_ab, ket_tensor, bra_det, ket_det );
  new_sigma_data_map_->emplace( "sigma_ab_vb", sigma_ab );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_ab_vb( shared_ptr<Vector_Bundle<DataType>> sigma_ab,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_ab_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  assert ( bra_det->nelea() == ket_det->nelea()+1  && bra_det->neleb() == ket_det->neleb()-1 );
  assert( sigma_ab->index_range_vec_[0].range().size() == 1 ) ; // TODO : catch for split CI block; must fix, but OK for now (for loops over blocks redundant)
  const int norb = ket_det->norb();
  const int ket_lenb = ket_det->lenb();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();


  SMITH::IndexRange ket_idx_range = ket_tensor->indexrange().front();
  for ( SMITH::Index block : ket_idx_range.range() )  {
    vector<SMITH::Index> ket_block_id = { block };
    unique_ptr<DataType[]> ket_data = ket_tensor->get_block( ket_block_id); 
    DataType* ket_ptr_orig = ket_data.get();
   
    // Slow, but simple to check and parallelize
    for (int ii = 0; ii != norb; ii++) {
      for (int jj = 0; jj != norb; jj++) {
    
        vector<int> orb_ids = { ii, jj };
        shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_ab->get_new_vector( true );

        bool ij_vec_sparse = true;
        for ( SMITH::Index& ij_block_id : sigma_ab->index_range_vec_[0] ) {
          
          DataType* ket_ptr = ket_ptr_orig;
          vector<SMITH::Index> sigma_ij_block_id = { ij_block_id };
          unique_ptr<DataType[]> sigma_ij_block_data = sigma_ij_vec->get_block( sigma_ij_block_id );
          DataType* sigma_ab_ij_ptr = sigma_ij_block_data.get();
          

          for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
            for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it, ++ket_ptr ) {
              if ( !(*abit_it)[ii] && (*bbit_it)[jj] ) {
          
                bitset<nbit__> bbit_bra = *bbit_it;
                DataType op_phase =  ket_det->sign<1>( bbit_bra, jj);
                bbit_bra.reset(jj);
          
                bitset<nbit__> abit_bra = *abit_it;
                op_phase *= ket_det->sign<0>( abit_bra, ii);
                abit_bra.set(ii);
                {
                  DataType* bra_ptr = sigma_ab_ij_ptr + bra_det->lexical<0>( abit_bra ) * bra_lenb + bra_det->lexical<1>( bbit_bra ); 
                  *bra_ptr += op_phase* (*ket_ptr);
                }
              }
            }
          }
          
          if ( ( abs(*(max_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ ) ||
               ( abs(*(min_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ )    ) { 
            ij_vec_sparse = false;
            sigma_ij_vec->put_block( sigma_ij_block_data, sigma_ij_block_id ); 
          }
        }
        sigma_ab->set_sparsity( orb_ids, ij_vec_sparse );
        if ( !ij_vec_sparse )
          sigma_ab->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
      }
    }
  } // end ket block loop
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ba_vb( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_ba_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = det_old_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = det_old_map_->at(bra_name);

  { // TEST
  bra_det = make_shared<Determinants>( ket_det->norb(), ket_det->nelea()-1, ket_det->neleb()+1, false /*compress*/, true /*mute*/);
  } // ENDTEST

  shared_ptr<Vector_Bundle<DataType>> sigma_ba;
  if ( new_sigma) {
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_ba = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );
  } else { 
    sigma_ba = new_sigma_data_map_->at( gamma_info->sigma_name() ) ;
  }

  convert_civec_to_tensor( ket_name ); 
  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name); 

  sigma2_ba_vb( sigma_ba, ket_tensor, bra_det, ket_det );
  new_sigma_data_map_->emplace( "sigma_ba_vb", sigma_ba );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_ba_vb( shared_ptr<Vector_Bundle<DataType>> sigma_ba,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_ba_vb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  assert ( bra_det->nelea() == ket_det->nelea()-1  && bra_det->neleb() == ket_det->neleb()+1 );
  const int norb = ket_det->norb();
  const int ket_lenb = ket_det->lenb();
  const int bra_lenb = bra_det->lenb();
  size_t bra_length = bra_det->lena()*bra_det->lenb();

  SMITH::IndexRange ket_idx_range = ket_tensor->indexrange().front();

  assert( sigma_ba->index_range_vec_[0].range().size() == 1 ) ; // TODO : catches for split CI block; must fix, but OK for now (for loops over blocks redundant)
  assert( ket_idx_range.range().size() == 1 ) ;
 
  for ( SMITH::Index block : ket_idx_range.range() )  {
    vector<SMITH::Index> ket_block_id = { block };
    unique_ptr<DataType[]> ket_data = ket_tensor->get_block( ket_block_id ); 
    DataType* ket_ptr_orig = ket_data.get();
   
    // Slow, but simple to check and parallelize
    for (int ii = 0; ii != norb; ii++) {
      for (int jj = 0; jj != norb; jj++) {
    
        vector<int> orb_ids = { ii, jj };
        shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_ba->get_new_vector( true );

        bool ij_vec_sparse = true;
        for ( SMITH::Index& ij_block_id : sigma_ba->index_range_vec_[0] ) {
          vector<SMITH::Index> sigma_ij_block_id = { ij_block_id };
          unique_ptr<DataType[]> sigma_ij_block_data = sigma_ij_vec->get_block( sigma_ij_block_id );
          DataType* sigma_ba_ij_ptr = sigma_ij_block_data.get();

          DataType* ket_ptr = ket_ptr_orig;
          for ( vector<bitset<nbit__>>::const_iterator abit_it = ket_det->string_bits_a().begin(); abit_it != ket_det->string_bits_a().end(); ++abit_it ) {
            for ( vector<bitset<nbit__>>::const_iterator bbit_it = ket_det->string_bits_b().begin(); bbit_it != ket_det->string_bits_b().end(); ++bbit_it, ++ket_ptr ) {
              if ( !(*bbit_it)[ii] && (*abit_it)[jj] ) {
          
                bitset<nbit__> abit_bra = *abit_it;
                DataType op_phase =  ket_det->sign<0>( abit_bra, jj);
                abit_bra.reset(jj);
          
                bitset<nbit__> bbit_bra = *bbit_it;
                op_phase *= ket_det->sign<1>( bbit_bra, ii);
                bbit_bra.set(ii);
                { 
                  DataType* bra_ptr = sigma_ba_ij_ptr + bra_det->lexical<0>( abit_bra ) * bra_lenb + bra_det->lexical<1>( bbit_bra ); 
                  *bra_ptr += op_phase* (*ket_ptr);
                }
              }
            }
          }
          
          if ( ( abs(*(max_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ ) ||
               ( abs(*(min_element(sigma_ij_block_data.get() , sigma_ij_block_data.get() + sigma_ij_block_id.front().size() ))) > thresh_ )    ) { 
            ij_vec_sparse = false;
            sigma_ij_vec->put_block( sigma_ij_block_data, sigma_ij_block_id ); 
          }
        }
        sigma_ba->set_sparsity( orb_ids, ij_vec_sparse );
        if ( !ij_vec_sparse )
          sigma_ba->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
      }
    }
  } // end ket block loop
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<vector<IndexRange>> B_Gamma_Computer::B_Gamma_Computer<DataType>::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(*ranges_str, "ranges_str" ) ; cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map_->at(rng));

  return ranges_Bagel;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma2_test( shared_ptr<GammaInfo_Base> gamma2_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::sigma2_test : " << gamma2_info->name() << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////
 
  string ket_name = gamma2_info->Ket_name();
  get_wfn_data( gamma2_info->Ket_info() );

  shared_ptr<Determinants> ket_det = det_old_map_->at( ket_name );
  shared_ptr<Civec>        ket = cvec_old_map_->at( ket_name );

  int norb = ket_det->norb();
  shared_ptr<Dvec> sigma2;

  if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha() ) {

    sigma2 = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );
    sigma_2a1( ket->data(), sigma2->data(0)->data(), ket_det );
    sigma_2a2( ket->data(), sigma2->data(0)->data(), ket_det );

    { //TEST

      {
      shared_ptr<Dvec> sigma2_aa_orig = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );
      sigma_2a1( ket->data(), sigma2_aa_orig->data(0)->data(), ket_det );

      sigma_aa_test( gamma2_info, true );
      shared_ptr<Dvec> sigma2_aa_test = dvec_sigma_map_->at("sigma_aa_test");
      vector<tuple<int,int,int>> failed_ids(0);
      cout << endl <<  "--------sigma2_aa_orig - sigma_test---------" << endl;
      for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
        for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
          for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
            DataType saa_diff = sigma2_aa_orig->data(rr*ket_det->norb()+ss)->data(qq) - sigma2_aa_test->data(rr*ket_det->norb()+ss)->data(qq) ;
            cout << saa_diff << "    "; cout.flush();
            if ( saa_diff != 0 ) {  
              int rrnss = rr*norb+ss;
              failed_ids.push_back( tie(rr, ss, rrnss )  );
            }
          }
        cout << endl;
      }

      {
      sigma_aa_vb( gamma2_info, true );
      shared_ptr<Vector_Bundle<DataType>> sigma_aa_vb_method =  new_sigma_data_map_->at("sigma_aa_vb_test"); 
      cout << endl <<  "--------sigma2_aa_orig - sigma_vb_test---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_aa_orig->data(rr + norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_aa_vb_method->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr - *ij_orig_ptr << "    ";
          sigma_aa_vb_method->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }

      }

      if ( failed_ids.size() > 0 ) {
        cout << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!  sigma2_aa FAILED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << endl <<  "--------sigma2_aa_orig ---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              DataType saa_diff = sigma2_aa_orig->data(rr*ket_det->norb()+ss)->data(qq);
              cout << saa_diff << "    "; cout.flush();
            }
          cout << endl;
        }
        cout << endl;
        cout << endl <<  "--------sigma_aa_test---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              DataType saa_diff = sigma2_aa_test->data(rr*ket_det->norb()+ss)->data(qq) ;
              cout << saa_diff << "    "; cout.flush();
            }
          cout << endl;
        }
        cout << endl << "++++++++++++++++++++ failed id pos ++++++++++++++++++++++++" << endl;
        for ( tuple<int,int,int>& fid : failed_ids ) 
          cout << "[ " << get<0>(fid) << " " << get<1>(fid) << " " << get<2>(fid) << " ]" << endl; 
 
      }
      }
      {
      convert_civec_to_tensor( ket_name );
      shared_ptr<Dvec> sigma2_bb_orig = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );
      sigma_2a2( ket->data(), sigma2_bb_orig->data(0)->data(), ket_det );

      sigma_bb_test( gamma2_info, true );
      shared_ptr<Dvec> sigma2_bb_test = dvec_sigma_map_->at("sigma_bb_test");
      vector<tuple<int,int,int>> failed_ids(0);
      cout << endl <<  "--------sigma2_bb_orig - sigma_bb_test---------" << endl;
      for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
        for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
          for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
            DataType sbb_diff = sigma2_bb_orig->data(rr*ket_det->norb()+ss)->data(qq) - sigma2_bb_test->data(rr*ket_det->norb()+ss)->data(qq) ;
            cout << sbb_diff << "    "; cout.flush();
            if ( sbb_diff != 0 ) { 
              int rrnss = rr*norb+ss;
              failed_ids.push_back( tie(rr, ss, rrnss )  );
            }
          }
        cout << endl;
      }
      {
      sigma_bb_vb( gamma2_info, true );
      shared_ptr<Vector_Bundle<DataType>> sigma_bb_vb =  new_sigma_data_map_->at("sigma_bb_vb_test"); 
      cout << endl <<  "--------sigma2_bb_orig - sigma_bb_vb---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_bb_orig->data(rr + norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_bb_vb->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr - *ij_orig_ptr << "    ";
          sigma_bb_vb->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }

      }


      if ( failed_ids.size() > 0 ) {
        cout << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!  sigma2_bb FAILED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        cout << endl <<  "--------sigma2_bb_orig ---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              DataType sbb_diff = sigma2_bb_orig->data(rr*ket_det->norb()+ss)->data(qq);
              cout << sbb_diff << "    "; cout.flush();
            }
          cout << endl;
        }
        cout << endl;
        cout << endl <<  "--------sigma_bb_test---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              DataType sbb_diff = sigma2_bb_test->data(rr*ket_det->norb()+ss)->data(qq) ;
              cout << sbb_diff << "    "; cout.flush();
            }
          cout << endl;
        }
        cout << endl << "++++++++++++++++++++ failed id pos ++++++++++++++++++++++++" << endl;
        for ( tuple<int,int,int>& fid : failed_ids ) 
          cout << "[ " << get<0>(fid) << " " << get<1>(fid) << " " << get<2>(fid) << " ]" << endl; 
 

      }
      }

      { 
        sigma_ab_test( gamma2_info, true ) ;
        shared_ptr<Dvec> sigma2_ab = dvec_sigma_map_->at("sigma_ab_test");
        cout << endl <<  "-------- sigma_ab_test---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              cout << sigma2_ab->data(rr*ket_det->norb()+ss)->data(qq) << "    "; cout.flush();
            }
          cout << endl;
        }
        sigma_ba_test( gamma2_info, true ) ;
        shared_ptr<Dvec> sigma2_ba = dvec_sigma_map_->at("sigma_ba_test");
        cout << endl <<  "-------- sigma_ba_test---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              cout << sigma2_ba->data(rr * ket_det->norb() + ss )->data(qq) << "    "; cout.flush();
            }
          cout << endl;
        }
        cout << endl <<  "-------- sigma_ba_test_ji + sigma_ab_test_ij ---------" << endl;
        for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
          for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
            for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
              cout << sigma2_ba->data(rr * ket_det->norb() + ss )->data(qq) - sigma2_ab->data(rr * ket_det->norb()+ ss)->data(qq) << "    "; cout.flush();
            }
          cout << endl;
        }
      }
    
      {
 
      shared_ptr<Dvec> sigma2_ab_orig = dvec_sigma_map_->at("sigma_ab_test");
      shared_ptr<const Determinants> bra_det = sigma2_ab_orig->det();      
      cout << endl <<  "--------sigma2_ab_orig ---------" << endl;
      for ( int rr = 0 ; rr != norb; rr++ ) 
        for ( int ss = 0 ; ss != norb; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ab_orig->data(rr+norb*ss)->data(); 
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_orig_ptr++ )
            cout << *ij_orig_ptr << "    ";
          cout << endl;
        }

      sigma_ab_vb( gamma2_info, true );
      shared_ptr<Vector_Bundle<DataType>> sigma_ab_vb =  new_sigma_data_map_->at("sigma_ab_vb"); 
      cout << endl <<  "-------- sigma_ab_vb---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ab_orig->data(rr+norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_ab_vb->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr << "    ";
          sigma_ab_vb->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }

    
      cout << endl <<  "-------- sigma_ab_vb - sigma_ab_orig---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ab_orig->data(rr+norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_ab_vb->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr - *ij_orig_ptr << "    ";
          sigma_ab_vb->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }
      }
      {
 
      shared_ptr<Dvec> sigma2_ba_orig = dvec_sigma_map_->at("sigma_ba_test");
      shared_ptr<const Determinants> bra_det = sigma2_ba_orig->det();      
      cout << endl <<  "--------sigma2_ba_orig ---------" << endl;
      for ( int rr = 0 ; rr != norb; rr++ ) 
        for ( int ss = 0 ; ss != norb; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ba_orig->data(rr+norb*ss)->data(); 
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_orig_ptr++ )
            cout << *ij_orig_ptr << "    ";
          cout << endl;
        }

      sigma_ba_vb( gamma2_info, true );
      shared_ptr<Vector_Bundle<DataType>> sigma_ba_vb =  new_sigma_data_map_->at("sigma_ba_vb"); 
      cout << endl <<  "-------- sigma_ba_vb---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ba_orig->data(rr+norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_ba_vb->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr << "    ";
          sigma_ba_vb->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }

    
      cout << endl <<  "-------- sigma_ba_vb - sigma_ba_orig---------" << endl;
      for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
        for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
          DataType* ij_orig_ptr = sigma2_ba_orig->data(rr+norb*ss)->data(); 
          std::unique_ptr<DataType[]> ij_vec = sigma_ba_vb->get_vector_ij_block( {rr, ss}, 0 );
          DataType* ij_vb_ptr = ij_vec.get();
          for ( int qq = 0 ; qq != bra_det->lena()*bra_det->lenb() ; qq++, ij_vb_ptr++, ij_orig_ptr++ )
            cout << *ij_vb_ptr - *ij_orig_ptr << "    ";
          sigma_ba_vb->put_vector_ij_block({rr,ss}, ij_vec, 0 );
          cout << endl;
        }
      }

    }
  }
  return; 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class B_Gamma_Computer::B_Gamma_Computer<double>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
