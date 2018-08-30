#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/gamma_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/util/prim_op.h>
#include <src/prop/proptool/debugging_utils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;
//#define __DEBUG_PROPTOOL_GAMMA_COMPUTER
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
Gamma_Computer::Gamma_Computer<DataType>::Gamma_Computer() {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::Gamma_Computer" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  civec_maxtile_ = 100000;// TODO set these elsewhere
  thresh_ = 1.0e-12;

  tensor_calc_           = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>();
  bagel_determinant_map_ = make_shared<std::map< std::string, std::shared_ptr<Determinants>>>();
  range_conversion_map_ =  make_shared<map< string, shared_ptr<SMITH::IndexRange>>>();
  gamma_info_map_ =  make_shared<map< string, shared_ptr<GammaInfo_Base>>>();
  civec_data_map_ =  make_shared<map< string, shared_ptr<Tensor_<DataType>>>>();

} 
/////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::get_gamma( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::get_gamma :  " <<  gamma_name << endl;
#endif //////////////////////////////////////////////////////////////////////////////

  if ( gamma_data_map_->find( gamma_name ) == gamma_data_map_->end()) {
    shared_ptr<GammaInfo_Base> gamma_info = gamma_info_map_->at(gamma_name);

    //note this has reverse iterators!
    if (gamma_info->order() > 2 ) { 
      compute_sigmaN( gamma_info );
    } else {
      compute_sigma2( gamma_info );
    }
    get_gammaN_from_sigmaN ( gamma_info ) ;
  }
  return;                              
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::compute_sigma2( shared_ptr<GammaInfo_Base> gamma2_info ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::compute_sigma2 : " << gamma2_info->name() << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  string ket_name = gamma2_info->Ket_name();

  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at( ket_name ); 

  if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha() ) { 
   sigma_aa( gamma2_info, true  );
   sigma_bb( gamma2_info, false );
    
  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()+1 ) { 
   sigma_ab( gamma2_info, true );

  } else if ( gamma2_info->Bra_nalpha() == gamma2_info->Ket_nalpha()-1 ) { 
   sigma_ba( gamma2_info, true );

  } else {
    throw logic_error( "this sigma: " + gamma2_info->sigma_name() + " is not implemented" ); 
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::get_gammaN_from_sigmaN( shared_ptr<GammaInfo_Base> gamma_n_info )  {  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::get_gammaN_from_sigmaN : " << gamma_n_info->name() << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////

  int order = gamma_n_info->order(); 
   
  shared_ptr<Vector_Bundle<DataType>> sigma_n = sigma_data_map_->at( gamma_n_info->sigma_name() ); 
  shared_ptr<Tensor_<DataType>> bra = civec_data_map_->at(gamma_n_info->Bra_name());
 
  vector<IndexRange> id_ranges = get_bagel_indexranges( *(gamma_n_info->id_ranges()) ); 
  shared_ptr<Tensor_<DataType>> gamma_n_tens = make_shared<Tensor_<DataType>>( id_ranges ); 
  gamma_n_tens->allocate();
  gamma_n_tens->zero();
  gamma_data_map_->emplace( gamma_n_info->name(), gamma_n_tens ); 
  
  vector<int> range_lengths = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(order,0);  
  vector<int> mins(order,0);  
  
  vector<vector<int>> block_offsets = get_block_offsets( id_ranges );
  
  do {
    vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
    vector<int> block_start = WickUtils::get_1d_from_2d( block_offsets, block_pos ); 
    vector<int> gamma_ids = block_start;
    vector<int> block_end = block_start;
      
    {
      vector<Index>::iterator ib_it = id_blocks.begin();
      for( vector<int>::iterator be_it = block_end.begin(); be_it != block_end.end(); ++be_it, ++ib_it ) 
        *be_it += ib_it->size()-1;
    }
    
    unique_ptr<DataType[]> gamma_block = gamma_n_tens->get_block( id_blocks );
    DataType* gamma_block_ptr = gamma_block.get(); 
  
    do {
      
     shared_ptr<Tensor_<DataType>> ket = sigma_n->vector_map(gamma_ids); 
     *gamma_block_ptr = ket->dot_product(bra);
     ++gamma_block_ptr;
  
    } while (fvec_cycle_skipper(gamma_ids, block_end, block_start ));
  
    gamma_n_tens->put_block(gamma_block, id_blocks);
  } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::compute_sigmaN( shared_ptr<GammaInfo_Base> gammaN_info )  {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::compute_sigmaN : " << gammaN_info->sigma_name() << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////

  int sorder  = gammaN_info->order(); 
  // check to see if previous sigma has been calculated, if not calculate it (recursive call here) 
  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at( gammaN_info->Bra_name() );  
  shared_ptr<Determinants> bra_det = bagel_determinant_map_->at( gammaN_info->Prev_Bra_name() );  

  size_t bra_length = bra_det->lena()*bra_det->lenb();
  int norb = ket_det->norb();
 
  vector<int> orb_ranges_sigma_n(sorder, norb );

  shared_ptr<Vector_Bundle<DataType>> sigma_n = make_shared<Vector_Bundle<DataType>>( orb_ranges_sigma_n, bra_length, civec_maxtile_, true, true, true );
 
  sigma_data_map_->emplace( gammaN_info->sigma_name(), sigma_n );

  // TODO should define norb so can be variable from gamma_info..
  vector<int> maxs_prev_sigma(sorder-2, norb-1);
  vector<int> mins_prev_sigma(sorder-2, 0);
  vector<int> orb_ids_prev_sigma = mins_prev_sigma;

  shared_ptr<Vector_Bundle<DataType>> prev_sigma = sigma_data_map_->at( gammaN_info->prev_sigma_name() );
  
  vector<bool> sigma_overwrite_pattern(sorder,false);
  sigma_overwrite_pattern[sorder-1] = true;
  sigma_overwrite_pattern[sorder-2] = true; 

  if ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha() ){ 

    do {  
      vector<int> maxs_sigma2(2, norb-1);
      shared_ptr<Tensor_<DataType>> ket = prev_sigma->vector_map(orb_ids_prev_sigma);
      auto tmp_sigma = make_shared<Vector_Bundle<DataType>>( maxs_sigma2, bra_length, civec_maxtile_, true, false, false );
      compute_eiej_on_ket( tmp_sigma, ket, bra_det, ket_det, "AA" ); 
      compute_eiej_on_ket( tmp_sigma, ket, bra_det, ket_det, "BB" ); 
      sigma_n->merge_fixed_ids( tmp_sigma, orb_ids_prev_sigma, sigma_overwrite_pattern, 'O' );
      
    } while(fvec_cycle_skipper( orb_ids_prev_sigma, maxs_prev_sigma, mins_prev_sigma) );
 
  } else if ( ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()+1 ) &&
              ( gammaN_info->Bra_nbeta() ==  gammaN_info->prev_Bra_nbeta()-1 ) ){ 
    do {  

      vector<int> maxs_sigma2(2, norb-1);
      vector<int> mins_sigma2(2, 0);
      shared_ptr<Tensor_<DataType>> ket = prev_sigma->vector_map(orb_ids_prev_sigma);

      auto tmp_sigma = make_shared<Vector_Bundle<DataType>>( maxs_sigma2, bra_length, civec_maxtile_, true, false, false );
      compute_eiej_on_ket( tmp_sigma, ket, bra_det, ket_det, "AB" ); 
      sigma_n-> merge_fixed_ids( tmp_sigma, orb_ids_prev_sigma, sigma_overwrite_pattern, 'O' );

    } while(fvec_cycle_skipper( orb_ids_prev_sigma, maxs_prev_sigma, mins_prev_sigma) );

  } else if ( ( gammaN_info->Bra_nalpha() ==  gammaN_info->prev_Bra_nalpha()-1 ) &&
              ( gammaN_info->Bra_nbeta()  ==  gammaN_info->prev_Bra_nbeta()+1 ) ){ 
   
    do {  
      vector<int> maxs_sigma2(2, norb-1);
      shared_ptr<Tensor_<DataType>> ket = prev_sigma->vector_map(orb_ids_prev_sigma);

      auto tmp_sigma = make_shared<Vector_Bundle<DataType>>( maxs_sigma2, bra_length, civec_maxtile_, true, false, false );
      compute_eiej_on_ket( tmp_sigma, ket, bra_det, ket_det, "BA" ); 
      sigma_n-> merge_fixed_ids( tmp_sigma, orb_ids_prev_sigma, sigma_overwrite_pattern, 'O' );

    } while(fvec_cycle_skipper( orb_ids_prev_sigma, maxs_prev_sigma, mins_prev_sigma) );

  } else { 
    throw logic_error( gammaN_info->prev_gamma_name() + " is not yet implemented! Aborting!!" );
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Gamma_Computer::Gamma_Computer<DataType>::compute_eiej_on_ket( shared_ptr<Vector_Bundle<DataType>> eiej_on_ket,
                                                                   shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                                   shared_ptr<Determinants> bra_det,
                                                                   shared_ptr<Determinants> ket_det,
                                                                   string transition_name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::compute_eiej_on_ket ";
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( transition_name == "AA") {
    assert( (bra_det->nelea() == ket_det->nelea()) && ( bra_det->neleb() == ket_det->neleb()));
    sigma2_aa( eiej_on_ket, ket_tensor, bra_det, ket_det );
  
  } else if ( transition_name == "BB") {
    assert( (bra_det->nelea() == ket_det->nelea()) && ( bra_det->neleb() == ket_det->neleb()));
    sigma2_bb( eiej_on_ket, ket_tensor, bra_det, ket_det );
    
  } else if ( transition_name == "AB") {
    assert( (bra_det->nelea()+1 == ket_det->nelea()) && ( bra_det->neleb()-1 == ket_det->neleb()));
    sigma2_ab( eiej_on_ket, ket_tensor, bra_det, ket_det );

  } else if ( transition_name == "BA") {
    assert( (bra_det->nelea()-1 == ket_det->nelea()) && ( bra_det->neleb()+1 == ket_det->neleb()));
    sigma2_ba( eiej_on_ket, ket_tensor, bra_det, ket_det );

  } else {
    throw logic_error( "Gamma_Computer::compute_eiej_on_ket : Aborted as this sigma is not implemented; not aa , bb, ab or ba" ); 
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::sigma_aa( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma_aa_test" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = bagel_determinant_map_->at(bra_name);

  shared_ptr<Vector_Bundle<DataType>> sigma_aa;
  auto sigma_map_loc = sigma_data_map_->find( gamma_info->sigma_name() );
  if ( new_sigma || sigma_map_loc == sigma_data_map_->end() ){ 
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_aa = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );
    if ( sigma_map_loc == sigma_data_map_->end() ){ 
      sigma_data_map_->insert( make_pair( gamma_info->sigma_name(), sigma_aa ) );
    } else { 
      sigma_map_loc->second = sigma_aa; 
    }
  } else {
    sigma_aa = sigma_map_loc->second;
  }
  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name);
  
  sigma2_aa( sigma_aa, ket_tensor, bra_det, ket_det );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Gamma_Computer::Gamma_Computer<DataType>::sigma2_aa( shared_ptr<Vector_Bundle<DataType>> sigma_aa,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma2_aa" << endl;
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
        
        if ( sigma_aa->sparsity_map(orb_ids) ) { 
          sigma_aa->set_sparsity( orb_ids, ij_vec_sparse ); 
          if ( !ij_vec_sparse )
            sigma_aa->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
        } else if (!ij_vec_sparse) { 
          sigma_aa->vector_map(orb_ids)->ax_plus_y(1.0, sigma_ij_vec); 
        }

      }
    }
    ket_tensor->put_block( ket_data );
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::sigma_bb( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma_bb" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = ket_det;

  shared_ptr<Vector_Bundle<DataType>> sigma_bb;
  auto sigma_map_loc = sigma_data_map_->find( gamma_info->sigma_name() );
  if ( new_sigma || sigma_map_loc == sigma_data_map_->end() ){ 

    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_bb = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );

    if ( sigma_map_loc == sigma_data_map_->end() ){ 
      sigma_data_map_->insert( make_pair( gamma_info->sigma_name(), sigma_bb ) );

    } else { 
      sigma_map_loc->second = sigma_bb; 

    }

  } else {
    sigma_bb = sigma_map_loc->second;

  }

  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name);
  
  sigma2_bb( sigma_bb, ket_tensor, bra_det, ket_det ) ;

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Gamma_Computer::Gamma_Computer<DataType>::sigma2_bb( shared_ptr<Vector_Bundle<DataType>> sigma_bb,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma2_bb" << endl;
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
          
          // TODO this is extremely inefficient, but all these sigma routines should be redone....
          unique_ptr<DataType[]> ij_block_t( new DataType[ij_block_id.size()] ); 
          fill_n( ij_block_t.get(), ij_block_id.size(), (DataType)(0.0) ); 
          vector<SMITH::Index> ij_block_id_vec = { ij_block_id };
      
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
            sigma_ij_vec->put_block( ij_block, ij_block_id_vec ); 

          }

        }
        
        if ( sigma_bb->sparsity_map(orb_ids) ) { 
          sigma_bb->set_sparsity( orb_ids, ij_vec_sparse ); 
          if ( !ij_vec_sparse )
            sigma_bb->set_vector( orb_ids , sigma_ij_vec, /*overwrite*/ true ); 
        } else if (!ij_vec_sparse) { 
          sigma_bb->vector_map(orb_ids)->ax_plus_y(1.0, sigma_ij_vec); 
        }
      }
    }
  }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void Gamma_Computer::Gamma_Computer<DataType>::sigma_ab( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma_ab" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = bagel_determinant_map_->at(bra_name);
  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name); 

  shared_ptr<Vector_Bundle<DataType>> sigma_ab;
  auto sigma_map_loc = sigma_data_map_->find( gamma_info->sigma_name() );

  if ( new_sigma || (sigma_map_loc == sigma_data_map_->end()) ) { 
    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_ab = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );

    if ( sigma_map_loc == sigma_data_map_->end() ){ 
      sigma_data_map_->insert( make_pair( gamma_info->sigma_name(), sigma_ab ) );
    } else { 
      sigma_map_loc->second = sigma_ab; 
    }

  } else {
    sigma_ab = sigma_map_loc->second;
  }
  
  sigma2_ab( sigma_ab, ket_tensor, bra_det, ket_det );
  sigma_data_map_->emplace( gamma_info->sigma_name(), sigma_ab );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Gamma_Computer::Gamma_Computer<DataType>::sigma2_ab( shared_ptr<Vector_Bundle<DataType>> sigma_ab,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma2_ab" << endl;
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
        shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_ab->get_vector( orb_ids, true );

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
void Gamma_Computer::Gamma_Computer<DataType>::sigma_ba( shared_ptr<GammaInfo_Base> gamma_info, bool new_sigma ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma_ba" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = gamma_info->Bra_name();
  string ket_name = gamma_info->Ket_name();
  
  shared_ptr<Determinants> ket_det = bagel_determinant_map_->at(ket_name);
  shared_ptr<Determinants> bra_det = bagel_determinant_map_->at(bra_name);

  shared_ptr<Vector_Bundle<DataType>> sigma_ba;
  auto sigma_map_loc = sigma_data_map_->find( gamma_info->sigma_name() );

  if ( new_sigma || sigma_map_loc == sigma_data_map_->end() ){ 

    vector<int> orb_id_ranges = { bra_det->norb(), bra_det->norb() };
    sigma_ba = make_shared<Vector_Bundle<DataType>>( orb_id_ranges, bra_det->lena()*bra_det->lenb(), civec_maxtile_, true, true, true );

    if ( sigma_map_loc == sigma_data_map_->end() ){ 
      sigma_data_map_->insert( make_pair( gamma_info->sigma_name(), sigma_ba ) );
    } else { 
      sigma_map_loc->second = sigma_ba; 
    }

  } else {
    sigma_ba = sigma_map_loc->second;
  }

  shared_ptr<SMITH::Tensor_<DataType>> ket_tensor = civec_data_map_->at(ket_name); 

  sigma2_ba( sigma_ba, ket_tensor, bra_det, ket_det );

  sigma_data_map_->emplace( gamma_info->sigma_name(), sigma_ba );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Gamma_Computer::Gamma_Computer<DataType>::sigma2_ba( shared_ptr<Vector_Bundle<DataType>> sigma_ba,
                                                            shared_ptr<SMITH::Tensor_<DataType>> ket_tensor,
                                                            shared_ptr<Determinants> bra_det, shared_ptr<Determinants> ket_det ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::Gamma_Computer<DataType>::sigma2_ba" << endl;
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
        shared_ptr<SMITH::Tensor_<DataType>> sigma_ij_vec = sigma_ba->get_vector( orb_ids, true );

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
vector<IndexRange> Gamma_Computer::Gamma_Computer<DataType>::get_bagel_indexranges(const vector<string>& ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_COMPUTER
cout << "Gamma_Computer::get_bagel_indexranges 1arg "; print_vector(ranges_str, "ranges_str" ) ; cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> ranges_bgl(ranges_str.size());
  vector<string>::const_iterator rs_it = ranges_str.begin();
  for ( vector<IndexRange>::iterator rb_it = ranges_bgl.begin(); rb_it != ranges_bgl.end(); ++rb_it, ++rs_it ) 
    *rb_it = *(range_conversion_map_->at(*rs_it) );

  return ranges_bgl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Gamma_Computer::Gamma_Computer<double>;
//template class Gamma_Computer::Gamma_Computer<std::complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
