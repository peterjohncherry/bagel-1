#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO this version is non-relativistic, must switch back to relativistic version with tensors at some point
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
B_Gamma_Computer::B_Gamma_Computer<DataType>::B_Gamma_Computer( std::shared_ptr<const Dvec > civectors ) : cc_(civectors){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::B_Gamma_Computer" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  dvec_sigma_map = make_shared<std::map< std::string, std::shared_ptr<Dvec>>>();
  det_old_map    = make_shared<std::map< std::string, std::shared_ptr<Determinants>>>();
  cvec_old_map   = make_shared<std::map< std::string, std::shared_ptr<Civec>>>();

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
      det_old_map->emplace( get_det_name( 'a', det->nelea(), 'A' , det->neleb(), det->norb() ), det ); 
      det = plus_det; 
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Temporary function; ideally B_Gamma_Computer should be constructed inside System_Computer, but the dependence on the Dvecs is preventing this
// So for the time being construct B_gamma_computer with the old maps in proptool, feed this to system computer, and then alter the Gamma_info_map
// as appropriate. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
B_Gamma_Computer::B_Gamma_Computer<DataType>::set_maps( std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                                              std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo_Base>>> Gamma_info_map_in,
                                              std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> gamma_data_map,
                                              std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> sigma_data_map,  
                                              std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<DataType>>>> civec_data_map ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << " B_Gamma_Computer::B_Gamma_Computer<DataType>::Set_maps" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  range_conversion_map = range_conversion_map_in;
  Gamma_info_map = Gamma_info_map_in;
  gamma_data_map_ = gamma_data_map;
  sigma_data_map_ = sigma_data_map;
  civec_data_map_ = civec_data_map;
 
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

    shared_ptr<GammaInfo_Base> gamma_info = Gamma_info_map->at(gamma_name);

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
  shared_ptr<Dvec> Dvec_sigma = dvec_sigma_map->at( sigma_name ) ;

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
  shared_ptr<Dvec> sigmaN = dvec_sigma_map->at(gammaN_info->sigma_name());
  string Bra_name = gammaN_info->Bra_name();
  shared_ptr<Civec> Bra = cvec_old_map->at( Bra_name );      
  shared_ptr<Determinants> Bra_det = det_old_map->at( Bra_name ) ; 

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
  shared_ptr<Determinants> ket_det = det_old_map->at( ket_name );  

  int sorder  = gammaN_info->order(); 
  int orb_dim = pow(ket_det->norb(), sorder-2);
  int orb2    = ket_det->norb()*ket_det->norb();

  if ( dvec_sigma_map->find(gammaN_info->prev_sigma_name()) == dvec_sigma_map->end() ) {
    if ( sorder > 4 ) {
      compute_sigmaN( Gamma_info_map->at(gammaN_info->prev_gamma_name()) ); 
    } else { 
      compute_sigma2( Gamma_info_map->at(gammaN_info->prev_gamma_name()) ); 
    }
  } 
  shared_ptr<Dvec> prev_sigma = dvec_sigma_map->at( gammaN_info->prev_sigma_name() ); 
  shared_ptr<Dvec> sigmaN     = make_shared<Dvec>( ket_det, orb_dim*orb2  );

  if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha() ){ 

    for ( int  ii = 0; ii != orb_dim; ii++) {
      sigma_2a1( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), ket_det );
      sigma_2a2( prev_sigma->data(ii)->data(), sigmaN->data(ii*orb2)->data(), ket_det );
    }

    dvec_sigma_map->emplace( gammaN_info->sigma_name(), sigmaN );

  } else if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()+1 ){ 

    shared_ptr<Determinants> bra_det = det_old_map->at( bra_name );  
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
  string Bra_name    =  gamma2_info->Bra_name();
 
  convert_civec_to_tensor( Bra_name );

  shared_ptr<Tensor_<DataType>> Tens_gamma2 = Tensor_Calc->contract_tensor_with_vector( sigma_data_map_->at( sigma2_name ), civec_data_map_->at( Bra_name ), 0 );
  gamma_data_map_->emplace( gamma2_info->name(), Tens_gamma2 ); 

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

  shared_ptr<Determinants> ket_det = det_old_map->at( ket_name ); 
  shared_ptr<Civec>        ket = cvec_old_map->at( ket_name );
  
  shared_ptr<Dvec> sigma2;
    
  if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha() ) { 

    sigma2 = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );

    sigma_2a1( ket->data(), sigma2->data(0)->data(), ket_det );
    { //TEST

      cout << endl <<  "--------sigma_2a1---------" << endl;
      for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
        for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
          for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
            cout << sigma2->data(rr*ket_det->norb()+ss)->data(qq) << "    "; cout.flush();
          }
        cout << endl;
      }
      
      shared_ptr<Dvec> sigma2_test = make_shared<Dvec>( ket_det, ket_det->norb()*ket_det->norb() );
      cout << endl << endl << "---------- sigma_2aa_alt ---------" << endl;
      sigma_2aa_alt( ket->data(), sigma2_test->data(0)->data(), ket_det, ket_det );
      cout << endl << endl << "---------- sigma_2aa_alt ---------" << endl;
      for ( int qq = 0 ; qq != ket_det->lena()*ket_det->lenb() ; qq++ ) {
        for ( int rr = 0 ; rr != ket_det->norb() ; rr++ ) 
          for ( int ss = 0 ; ss != ket_det->norb() ; ss++ ) {
            cout << sigma2_test->data(rr*ket_det->norb()+ss)->data(qq) << "    "; cout.flush();
          }
        cout << endl;
      }


    }

    cout << endl;
    sigma_2a2( ket->data(), sigma2->data(0)->data(), ket_det );
    
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()+1 ) { 

     shared_ptr<Determinants> bra_det = det_old_map->at( gamma2_info->Bra_name() );

     sigma2 = make_shared<Dvec>( bra_det, bra_det->norb()*bra_det->norb() );

     sigma_2ab( ket->data(), sigma2->data(0)->data(), bra_det, ket_det );
  
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()-1 ) { 

    throw logic_error( "(-a,+b) sigmas not implemented yet. Aborting !!  ");  

  } else {

    cout << "this sigma : "; cout.flush(); cout << gamma2_info->sigma_name() << " is not implemented" << endl;
    throw logic_error( "Aborting !!  ");  

  }

  dvec_sigma_map->emplace( gamma2_info->sigma_name(), sigma2 );

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
void B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2aa_alt( DataType* cvec_ptr, DataType* sigma_ptr, 
                                                                  shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::B_Gamma_Computer<DataType>::sigma_2aa_alt" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //get pointer in target tensor
  size_t bra_civec_length = bra_det->lena()* bra_det->lenb();
  size_t bra_lenb = bra_det->lenb();
  size_t ket_lenb = ket_det->lenb();

  size_t norb = bra_det->norb();
  for (unsigned int ii =0 ; ii!=norb ; ii++) { 
    for (unsigned int jj =0 ; jj!=norb ; jj++) {
      DataType* sigma_ij_0_ptr = sigma_ptr + (ii*norb+jj)*bra_civec_length;

      cout << " ii = " << ii ; cout.flush() ; cout << "  jj = " << jj ; cout.flush(); cout << endl;
      cout << "*sigma_ij_0_ptr = "; cout.flush();  cout << *sigma_ij_0_ptr << endl;

      //get position in bra vec by doing the opposite of the operation you are going to do on the ket
      for ( auto& a1 : bra_det->phidowna(ii)){
        cout << "X1a" << endl;
        for ( auto& b1 : bra_det->phiupa(jj)  ){
          cout << "X1b" << endl;
          DataType* sigma_aj_bi_ptr = sigma_ij_0_ptr + (a1.source*ket_lenb + b1.source);
          cout << " *sigma_aj_ai_ptr  = " ; cout.flush(); cout << *sigma_aj_bi_ptr << endl;
          //get position in ket by standard approach
          for ( auto& a2 : ket_det->phiupa(ii)){     
            cout << "X2" << endl;
            for ( auto& b2 : ket_det->phidowna(jj)) { 
             cout << "X3" << endl;
              *sigma_aj_bi_ptr += *(cvec_ptr + a2.source*ket_lenb + b2.source);
             cout << "X4" << endl;
            }
          }
        }
      }
    }
  }
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
  auto cvec_det_loc = det_old_map->find(cvec_name); 
  if( cvec_det_loc == det_old_map->end()){
    int II = cvec_info->state_num();
    shared_ptr<const Determinants> det_cvec_orig = cc_->data(II)->det();
    shared_ptr<Determinants> det_cvec = make_shared<Determinants>(det_cvec_orig->norb(), det_cvec_orig->nelea(), det_cvec_orig->neleb(), false, /*mute=*/true);
    det_old_map->emplace( cvec_name, det_cvec);
    cvec_old_map->emplace( cvec_name, make_shared<Civec>(*(cc_->data(II))) );
  }  
  cout << "got " << cvec_name << " in cvec_old_map " <<endl; 
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
 
    vector<IndexRange> civec_idxrng(1, *(range_conversion_map->at(civec_name)) );
    cout <<" civec_name = " << civec_name << endl;
    shared_ptr<Tensor_<DataType>> civec_tensor = make_shared<Tensor_<DataType>>( civec_idxrng );
    civec_tensor->allocate();
    civec_tensor->zero();
    size_t idx_position = 0;
    cout << "looking for old civec " << civec_name << " ....." ; cout.flush();
    shared_ptr<Civec> civector = cvec_old_map->at(civec_name);
    cout << "found it!! " << endl;
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr<vector<IndexRange>> B_Gamma_Computer::B_Gamma_Computer<DataType>::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_B_GAMMA_COMPUTER
cout << "B_Gamma_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(*ranges_str, "ranges_str" ) ; cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class B_Gamma_Computer::B_Gamma_Computer<double>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
