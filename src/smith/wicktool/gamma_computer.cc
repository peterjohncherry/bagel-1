#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/gamma_computer.h>
#include <src/smith/wicktool/WickUtils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace Tensor_Arithmetic;
using namespace Tensor_Arithmetic_Utils;
using namespace WickUtils;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Gamma_Computer::Gamma_Computer::Gamma_Computer( shared_ptr< map< string, shared_ptr<GammaInfo>>>          Gamma_info_map_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    CIvec_data_map_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    Sigma_data_map_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    Gamma_data_map_in,
                                                shared_ptr< map< string, shared_ptr<const Determinants>>> Determinants_map_in,
						shared_ptr< map< string, shared_ptr<IndexRange>>>         range_conversion_map_in ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::Gamma_Computer" << endl;
  maxtile  = 10000;

  Gamma_info_map             = Gamma_info_map_in;             cout << "set Gamma_info_map            "<< endl; 
  CIvec_data_map       = CIvec_data_map_in;       cout << "set CIvec_data_map      "<< endl;  
  Sigma_data_map       = Sigma_data_map_in;       cout << "set Sigma_data_map      "<< endl;  
  Gamma_data_map       = Gamma_data_map_in;       cout << "set Gamma_data_map      "<< endl;  
  Determinants_map     = Determinants_map_in;     cout << "set Determinants_map    "<< endl; 
  range_conversion_map = range_conversion_map_in; cout << "set range_conversion_map"<< endl; 

  Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>();

  //tester();
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::get_gamma_tensor_test( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;


  if( Gamma_data_map->find(gamma_name) != Gamma_data_map->end()){ 

    cout << "already have data for " << gamma_name << endl;  
 
  } else { 
    
    shared_ptr<GammaInfo> gamma_info = Gamma_info_map->at(gamma_name) ;
    int order = gamma_info->order;
    //for now just use specialized routines, this must be made generic at some point
    if ( order == 2 ) { 

      build_gamma_2idx_tensor( gamma_name ) ;
      assert( gamma_2idx_contract_test( gamma_name ) );
      
    } else {
    
      build_gammaN_tensor( gamma_info ) ;

    }
  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_gammaN_tensor(shared_ptr<GammaInfo> gamma_info )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Gamma_Computer::build_gammaN_tensor" << endl;
   
   build_sigmaN_tensor(gamma_info);

   //TODO contract sigma 
   
   return;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigmaN_tensor(shared_ptr<GammaInfo> sigmaN_info )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// sn  :  sigmaN   ,  ps : prev_sigma    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigmaN_tensor" << endl;
 
  string sigma_name      = "S_"+sigmaN_info->name;
  string prev_sigma_name = "S_"+sigmaN_info->predecessor_gamma_name;  

  shared_ptr<CIVecInfo<double>> Bra_info = sigmaN_info->Bra_info;  
  shared_ptr<Tensor_<double>>   Bra = CIvec_data_map->at( Bra_info->name() );  

  shared_ptr<CIVecInfo<double>>  Ket_info   = sigmaN_info->Ket_info;  
  shared_ptr<Tensor_<double>>    prev_sigma = Sigma_data_map->at( prev_sigma_name );
  shared_ptr<const Determinants> rhs_det    = Determinants_map->at( Ket_info->name() ); 

  shared_ptr<GammaInfo> prev_sigma_info = Gamma_info_map->at( prev_sigma_name );

  int order = sigmaN_info->order;  

  vector<IndexRange> ranges_sn(order+1);
  for ( int ii = 0 ; ii != order ; ii++) 
    ranges_sn[ii] = *(range_conversion_map->at(sigmaN_info->id_ranges->at(ii)));
  ranges_sn.back() = *(range_conversion_map->at(Bra_info->name()));//fix this; ok in nr, but in rel lhs det space for gamma is not same as lhs det space for sigma

  shared_ptr<Tensor_<double>> sigmaN = make_shared<Tensor_<double>>( ranges_sn );

  vector<IndexRange> ranges_ps(order-1);
  for ( int ii = 0 ; ii != order-2 ; ii++) 
    ranges_ps[ii] = *(range_conversion_map->at(prev_sigma_info->id_ranges->at(ii)));
  ranges_ps.back() = *(range_conversion_map->at(Ket_info->name()));//fix this; ok in nr, but in rel lhs det space for gamma is not same as lhs det space for sigma



  shared_ptr<vector<int>> mins_ps          = make_shared<vector<int>>( ranges_ps.size(), 0 );  
  shared_ptr<vector<int>> block_pos_ps     = make_shared<vector<int>>( ranges_ps.size(), 0 );  
  shared_ptr<vector<int>> range_lengths_ps = get_range_lengths(ranges_ps); 

  shared_ptr<vector<vector<int>>> block_offsets_ps = get_block_offsets( ranges_ps ) ;
  vector<int> offsets_ps(order+1); 
  do {

    for ( int ii = 0 ; ii != offsets_ps.size(); ii++ )
      offsets_ps[ii] = block_offsets_ps->at(ii)[block_pos_ps->at(ii)]; 

    vector<Index> id_blocks_ps = *(get_rng_blocks( block_pos_ps, ranges_ps));

    vector<IndexRange> ranges_Iij = { ranges_sn[0], ranges_sn[1], ranges_sn.back(); }; 

    shared_ptr<vector<int>> mins_Iij          = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> block_pos_Iij     = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> range_lengths_Iij = get_range_lengths(ranges_Iij); 

    shared_ptr<vector<vector<int>>> block_offsets_Iij = get_block_offsets( ranges_Iij ) ;

    vector<int> offsets_Iij(2); 
    do { 

      vector<Index> id_blocks_Iij = *(get_rng_blocks( block_pos_Iij, ranges_Iij));
      for ( int ii = 0 ; ii != 3 ; ii++ ) 
        offsets_Iij[ii] = block_offsets_Iij->at(ii)[block_pos_Iij->at(ii)];    
   
      
    build_sigmaN_block( sigmaN, id_blocks_ps, offsets_ps, prev_sigma, id_blocks_Iij, offsets_Iij ) ;


    } while (fvec_cycle(block_pos_Iij, range_lengths_Iij, mins_Iij ));


  
  } while (fvec_cycle(block_pos_ps, range_lengths_ps, mins_ps ));

//  for ( int  ii = 0; ii != orb_dim; ii++) {

    // unique_ptr<double[]> get_sigma_part(predecessor_name, gamma_name); 
    // contract block with Bra to get gamma 

//  } 
 

  return;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_sigmaN_block( string sigmaN_name,     vector<Index>& id_blocks_ps, vector<int>& offsets_ps,
                                                    string prev_sigma_name, vector<Index>& id_blocks_Iij, vector<int>& offsets_Iij ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (Bra_name == Ket_name) { 

    shared_ptr<Tensor_<double>> sigmaN      = Sigma_data_map->at(sigmaN_name);
    shared_ptr<Tensor_<double>> prev_sigma  = Sigma_data_map->at(prev_sigma_name);

    unique_ptr<double[]> prev_sigma_block = prev_sigma->get_block( id_block_ps );
    
    vector<int> Iij_block_size      = id_blocks_Iij[0].size() * id_blocks_Iij[1].size()*id_blocks_Iij[2].size();
    vector<int> rhs_ci_block_size   = id_blocks_ps.back().size();
    vector<int> ps_orb_block_size = 1;
    for ( int bl_size : id_blocks_ps ) 
       ps_orb_block_size *= bl_size;

    unique_ptr<double[]> sigmaN_block( new double[ij_block_size*ps_orb_block_size*rhs_ci_block_size] );

    shared_ptr<Determinants> Ket_det  = Determinants_map->at(Sigma_data_map->at(sigmaN_name)->Ket_info->name());

    for ( int  ii = 0; ii != ps_orb_block_size; ii++) {
      sigma_2a1( prev_sigma_block.get()+(ii*rhs_ci_block_size), sigmaN_block.get()+(ii*Iij_block_size), Ket_det );
      sigma_2a2( prev_sigma_block.get()+(ii*rhs_ci_block_size), sigmaN_block.get()+(ii*Iij_block_size), Ket_det );
    }

  } else {

    cout << "spin flips not implemented yet " <<endl;
  }


}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::sigma_2a1(double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets,
                                               int cvec_offset , int sigma_cvec_offset, int ii_offset, int jj_offset ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                               
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
                                                                                                                          
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                             
void Gamma_Computer::Gamma_Computer::sigma_2a2( double* cvec_ptr, double* sigma_ptr, shared_ptr<Determinants> dets, 
                                                int cvec_offset , int sigma_cvec_offset, int ii_offset, int jj_offset ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////                             
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets the gammas in tensor format. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::get_gamma_tensor( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  if( Gamma_data_map->find(gamma_name) != Gamma_data_map->end()){ 

    cout << "already have data for " << gamma_name << endl;  
 
  } else { 
    
    //for now just use specialized routines, this must be made generic at some point
    if (Gamma_info_map->at(gamma_name)->id_ranges->size() == 2 ) { 

      build_gamma_2idx_tensor( gamma_name ) ;
      assert( gamma_2idx_contract_test( gamma_name ) );
      
    } else if (Gamma_info_map->at(gamma_name)->id_ranges->size() == 4 ) { 

      build_gamma_4idx_tensor( gamma_name );
      assert (gamma_4idx_contract_test( gamma_name ) );

      build_sigma_4idx_tensor( Gamma_info_map->at(gamma_name) );

    } else if (Gamma_info_map->at(gamma_name)->id_ranges->size() == 6 ) { 
       cout << " 6-index stuff not implemented, cannot calculate " << gamma_name << endl;

    } else if (Gamma_info_map->at(gamma_name)->id_ranges->size() == 8 ) { 
       cout << " 8-index stuff not implemented, cannot calculate " << gamma_name << endl;
    }    
  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_gamma_2idx_tensor( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "build_gamma_2idx_tensor : " << gamma_name << endl;

  shared_ptr<GammaInfo> gamma_info      =  Gamma_info_map->at(gamma_name);

  shared_ptr<CIVecInfo<double>> BraInfo =  gamma_info->Bra_info;
  shared_ptr<CIVecInfo<double>> KetInfo =  gamma_info->Ket_info;
  string Bra_name = BraInfo->name();  
  string Ket_name = KetInfo->name();  

  //norm check
  assert ( abs (1.0 - Tensor_Calc->contract_vectors(  CIvec_data_map->at(Bra_name) , CIvec_data_map->at(Bra_name)) ) < 0.0000001 ); 
 
  string sigma_name = "S_"+gamma_name;

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
     
     cout << "already got " << sigma_name << ", so use it " << endl;
  
  } else { //should replace this with blockwise and immediate contraction build; avoid storing sigma unnecessarily. 

    build_sigma_2idx_tensor( gamma_info );
  
    cout << "got sigma tensor :  " <<  sigma_name << endl;  

    shared_ptr<Tensor_<double>> gamma_2idx =  Tensor_Calc->contract_tensor_with_vector( Sigma_data_map->at(sigma_name), CIvec_data_map->at(Bra_name),  make_pair(2,0));
 
    Gamma_data_map->emplace( gamma_name, gamma_2idx );

  }
 
  return;
 
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_sigma_2idx_tensor(shared_ptr<GammaInfo> gamma_info )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_sigma_2idx_tensor" << endl;
 
  string sigma_name = "S_"+gamma_info->name;

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

  } else { 
    
    string Bra_name = gamma_info->Bra_info->name();
    string Ket_name = gamma_info->Ket_info->name();


    IndexRange Bra_range = CIvec_data_map->at(Bra_name)->indexrange()[0];
    
    shared_ptr<vector<IndexRange>> orb_ranges = Get_Bagel_IndexRanges(gamma_info->id_ranges );   
    
    auto sigma_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { orb_ranges->at(0), orb_ranges->at(1), Bra_range } );

    shared_ptr<Tensor_<double>> sigma_tensor = make_shared<Tensor_<double>>(*(sigma_ranges));
    sigma_tensor->allocate();
    sigma_tensor->zero();
   
    shared_ptr<vector<int>> mins          = make_shared<vector<int>>( sigma_ranges->size(), 0 );  
    shared_ptr<vector<int>> block_pos     = make_shared<vector<int>>( sigma_ranges->size(), 0 );  
    shared_ptr<vector<int>> range_lengths = get_range_lengths(sigma_ranges); 

    vector<int> sigma_offsets = { 0, 0, 0 };
    shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets( sigma_ranges ) ;

    // loop through sigma ranges,  loop over Ket ranges inside build_sigma_block;  
    do {

      for ( int ii = 0 ; ii != sigma_offsets.size(); ii++ )
        sigma_offsets[ii] = block_offsets->at(ii)[block_pos->at(ii)]; 

      vector<Index> sigma_id_blocks = *(get_rng_blocks( block_pos, *sigma_ranges));

      build_sigma_block( sigma_tensor, sigma_id_blocks, sigma_offsets, Ket_name ) ;
    
    } while (fvec_cycle(block_pos, range_lengths, mins ));
    
    Sigma_data_map->emplace(sigma_name, sigma_tensor); 

  }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigma_block( shared_ptr<Tensor_<double>> sigma_tensor, vector<Index>& sigma_id_blocks,
                                                        vector<int>& sigma_offsets, string Ket_name  ) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigma_block" << endl; 

  size_t iblock_size    = sigma_id_blocks[0].size();
  size_t jblock_size    = sigma_id_blocks[1].size();
  size_t Bra_block_size = sigma_id_blocks[2].size();
  const size_t sigma_block_size = Bra_block_size*iblock_size*jblock_size;

  shared_ptr<IndexRange>      Ket_range = range_conversion_map->at(Ket_name);
  shared_ptr<Tensor_<double>> Ket_Tens  = CIvec_data_map->at(Ket_name);

  unique_ptr<double[]> sigma_block(new double[sigma_block_size])  ;
  std::fill_n(sigma_block.get(), sigma_block_size, 0.0);
  double* sigma_ptr =sigma_block.get();

  shared_ptr<const Determinants> Ket_det = Determinants_map->at(Ket_name);
  const int lena = Ket_det->lena();
  const int lenb = Ket_det->lenb();

  size_t Ket_offset = 0;

  // Must be changed to use BLAS, however, you will have to be careful about the ranges; 
  // ci_block needs to have a length which is some integer multiple of lenb.
  for ( Index Ket_idx_block : Ket_range->range()) {  
    unique_ptr<double[]>  Ket_block = Ket_Tens->get_block( vector<Index> {Ket_idx_block} );
    
    double* Ket_ptr = Ket_block.get();
    size_t Ket_block_size = Ket_idx_block.size();
    
    for( size_t ii = sigma_offsets[0]; ii != iblock_size + sigma_offsets[0]; ii++ ) {
      for( size_t jj = sigma_offsets[1]; jj != jblock_size + sigma_offsets[1]; jj++ ) {
    
        for ( DetMap iter : Ket_det->phia(ii, jj)) {
   
          int Ket_shift = iter.target*lenb - Ket_offset; 
          if ( (Ket_shift < 0 ) || (Ket_shift >= Ket_block_size) )  continue;
    
          int Bra_shift = iter.source*lenb - sigma_offsets[2];
          if ( (Bra_shift < 0 ) || (Bra_shift >= Bra_block_size) )  continue;
             
          size_t loop_limit = lenb > Bra_block_size ? Bra_block_size : lenb;      //fix
          loop_limit = Ket_block_size > loop_limit ? loop_limit : Ket_block_size; //fix
          double sign = static_cast<double>(iter.sign);
          for( size_t ib = 0; ib != loop_limit; ib++) 
           *(sigma_ptr+Bra_shift+ib) += (*(Ket_ptr+Ket_shift+ib) * sign);
        }
    
        for( size_t ia = 0; ia != lena; ia++) {
          for ( DetMap iter : Ket_det->phib(ii, jj)) {
        
            int Ket_shift = iter.target + (ia*lenb) - Ket_offset;
            if ( (Ket_shift < 0 ) || (Ket_shift >= Ket_block_size) ) continue;
   
            int Bra_shift = iter.source + (ia*lenb) - sigma_offsets[2];
            if ( (Bra_shift < 0 ) || (Bra_shift >= Bra_block_size) ) continue;
            
            double sign = static_cast<double>(iter.sign);
        
            *(sigma_ptr + Bra_shift) += (*(Ket_ptr+Ket_shift) * sign);
          
          }
        } 
        sigma_ptr += Bra_block_size;
      }   
    }     
    Ket_offset += Ket_idx_block.size(); 
  }
  cout << " got a sigma block .... " ; cout.flush();
  sigma_tensor->add_block(sigma_block, sigma_id_blocks ) ; 
  cout << " and put it in the tensor! " << endl;
  return;
 
}       


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigma_4idx_tensor(shared_ptr<GammaInfo> gamma4_info )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_sigma4_tensor" << endl;

  string sigma4_name = "S_"+gamma4_info->name;

  if ( Sigma_data_map->find(sigma4_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma4_name << endl;
  
  } else {  

//    assert ( build_sigma4_tensor_tests( gamma4_info ) );

    //acquire sigma2_KJ
    shared_ptr<GammaInfo> sigma2_KJ_info = Gamma_info_map->at( gamma4_info->sub_gammas(1));
    build_sigma_2idx_tensor( sigma2_KJ_info );  

    shared_ptr<Tensor_<double>>    sigma2_KJ        = Sigma_data_map->at( "S_"+ sigma2_KJ_info->name ); 
    shared_ptr<vector<IndexRange>> sigma2_KJ_ranges = make_shared<vector<IndexRange>>( sigma2_KJ->indexrange() ); 

    //get ranges for sigma4_IJ and build    
    shared_ptr<vector<IndexRange>> orb_ranges  = Get_Bagel_IndexRanges( gamma4_info->id_ranges );   
    IndexRange IBra_idxrng = CIvec_data_map->at( sigma2_KJ_info->Bra_info->name() )->indexrange()[0]; 
    auto sigma4_IJ_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { orb_ranges->at(0), orb_ranges->at(1), orb_ranges->at(2), orb_ranges->at(3), IBra_idxrng} );

    shared_ptr<Tensor_<double>> sigma4_IJ = make_shared<Tensor_<double>>(*sigma4_IJ_ranges);
    sigma4_IJ->allocate();
    sigma4_IJ->zero();


    shared_ptr<vector<IndexRange>> ijI_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { orb_ranges->at(0), orb_ranges->at(1), IBra_idxrng } ); 
    shared_ptr<vector<int>> ijI_mins          = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> ijI_block_pos     = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> ijI_range_lengths = get_range_lengths(ijI_ranges); 
  
    shared_ptr<vector<int>> sigma2_KJ_mins          = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> sigma2_KJ_block_pos     = make_shared<vector<int>>( 3, 0 );  
    shared_ptr<vector<int>> sigma2_KJ_range_lengths = get_range_lengths(sigma2_KJ_ranges); 

    shared_ptr<vector<vector<int>>> ijI_block_offsets = get_block_offsets( ijI_ranges ) ;
    shared_ptr<vector<vector<int>>> sigma2_KJ_block_offsets = get_block_offsets( sigma2_KJ_ranges ) ;

    // really want to loop over the 4-electron sigma,  however, it is easier to have two seperate loops,
    // one for each two electron sigma, and then combine the blocks and offsets into the arguments we need.
    do {

      print_vector(*ijI_block_pos, "ijI_block_pos"); cout << endl;

      vector<int> ijI_offsets = { ijI_block_offsets->at(0)[ijI_block_pos->at(0)],  ijI_block_offsets->at(1)[ijI_block_pos->at(1)],
                                  ijI_block_offsets->at(2)[ijI_block_pos->at(2)] }; 

      vector<Index> ijI_id_blocks = *(get_rng_blocks( ijI_block_pos, *ijI_ranges));
 
      do {
        print_vector(*sigma2_KJ_block_pos, "sigma2_KJ_block_pos"); cout << endl;

        vector<int> sigma2_KJ_offsets(3);
        for ( int ii = 0 ; ii != sigma2_KJ_block_offsets->size(); ii++ )  
           sigma2_KJ_offsets[ii] = sigma2_KJ_block_offsets->at(ii)[sigma2_KJ_block_pos->at(ii)]; 
        
        vector<int> sigma4_IJ_offsets = { ijI_offsets[0], ijI_offsets[1], sigma2_KJ_offsets[0], sigma2_KJ_offsets[1], ijI_offsets[2]};

        print_vector( sigma4_IJ_offsets, "sigma4_IJ_offsets" ) ;

        vector<Index> sigma2_KJ_id_blocks = *(get_rng_blocks( sigma2_KJ_block_pos, *sigma2_KJ_ranges));

        shared_ptr<vector<int>> sigma4_IJ_block_pos =
                                make_shared<vector<int>>( vector<int> { ijI_block_pos->at(0), ijI_block_pos->at(1), sigma2_KJ_block_pos->at(0), sigma2_KJ_block_pos->at(1), ijI_block_pos->at(2)} );

        vector<Index> sigma4_IJ_id_blocks = *( get_rng_blocks( sigma4_IJ_block_pos, *sigma4_IJ_ranges ) ); 
 
        build_sigma_4idx_block_from_sigma_2idx_block( sigma4_IJ,  sigma4_IJ_id_blocks, sigma4_IJ_offsets,
                                                      sigma2_KJ,  sigma2_KJ_id_blocks, sigma2_KJ_offsets,
                                                      sigma2_KJ_info->Bra_info->name() )  ;

        
      } while (fvec_cycle_skipper(sigma2_KJ_block_pos, sigma2_KJ_range_lengths, sigma2_KJ_mins ));

    } while (fvec_cycle_skipper( ijI_block_pos, ijI_range_lengths, ijI_mins ));

    cout << endl << "----------------------------------------------------------------------------------------" << endl;
    Print_Tensor(sigma4_IJ, "sigma4_IJ"); 
    cout << endl << "----------------------------------------------------------------------------------------" << endl;
 
    shared_ptr<Tensor_<double>> gamma4_IJ =  Tensor_Calc->contract_tensor_with_vector( sigma4_IJ, CIvec_data_map->at(gamma4_info->Bra_info->name()),  make_pair(4,0));

    cout << endl << "----------------------------------------------------------------------------------------" << endl;
    Print_Tensor(gamma4_IJ, "gamma4_IJ from sigma4_IJ"); 
    cout << endl << "----------------------------------------------------------------------------------------" << endl;

    vector<int> all_ids = { 0,1,2,3} ;
    shared_ptr<Tensor_<double>> gamma4_IJ_full_contract = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor_new( gamma4_IJ, all_ids) ;
    cout << " gamma4_IJ_full_contract->rms() = " << gamma4_IJ_full_contract->rms() <<  endl;
 
    Sigma_data_map->emplace(sigma4_name, sigma4_IJ); 

  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// gets block of            sigma_{ijkl}^{I} =   \sum_{JK} < IBra | i^{+} j | KKet > < KBra |   k^{+} l | JKet > c_{J}
// this routine does   :   \sum_{K}  < IBra | i^{+} j | KKet > sigma_{kl}^{K}
// where sigma_{kl}^{K} =  \sum_{kl}\sum_{J} < KBra |   k^{+} l | JKet >
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_sigma_4idx_block_from_sigma_2idx_block( shared_ptr<Tensor_<double>> sigma4_IJ, vector<Index>& sigma4_IJ_id_blocks, vector<int>& sigma4_IJ_offsets,
                                                                              shared_ptr<Tensor_<double>> sigma2_KJ, vector<Index>& sigma2_KJ_id_blocks, vector<int>& sigma2_KJ_offsets,
                                                                              string KKet_name  )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigma4_block_from_sigma_2idx_block" << endl;

  size_t iblock_size     = sigma4_IJ_id_blocks[0].size();cout << " iblock_size     = "    << iblock_size     ; 
  size_t siblock_size    = sigma2_KJ_id_blocks[0].size();cout << " siblock_size     = " << siblock_size    << endl; 
  size_t jblock_size     = sigma4_IJ_id_blocks[1].size();cout << " jblock_size     = "    << jblock_size     ;
  size_t sjblock_size    = sigma2_KJ_id_blocks[1].size();cout << " sjblock_size     = " << sjblock_size    << endl;
  size_t kblock_size     = sigma4_IJ_id_blocks[2].size();cout << " kblock_size     = "    << kblock_size     << endl;
  size_t lblock_size     = sigma4_IJ_id_blocks[3].size();cout << " lblock_size     = "    << lblock_size     << endl;
                                                                                                            
  size_t IBra_block_size = sigma4_IJ_id_blocks[4].size();cout << " IBra_block_size = " << IBra_block_size << endl;
  size_t KKet_block_size = sigma2_KJ_id_blocks[2].size();cout << " KKet_block_size = " << KKet_block_size << endl;

  const size_t sigma4_IJ_block_size = IBra_block_size*iblock_size*jblock_size*kblock_size*lblock_size;

  shared_ptr<Tensor_<double>> gamma_2idx_IJ_test =  Tensor_Calc->contract_tensor_with_vector( sigma2_KJ, CIvec_data_map->at(KKet_name),  make_pair(2,0));
  Print_Tensor(gamma_2idx_IJ_test, "gamma_2idx_IJ_test");


  unique_ptr<double[]> sigma4_IJ_block(new double[sigma4_IJ_block_size])  ;
  std::fill_n(sigma4_IJ_block.get(), sigma4_IJ_block_size, 0.0);

  shared_ptr<const Determinants> KKet_det = Determinants_map->at( KKet_name );
  const int lena = KKet_det->lena();
  const int lenb = KKet_det->lenb();

  Index KKet_idx_block = sigma2_KJ_id_blocks[2] ;
  size_t KKet_offset   = sigma2_KJ_offsets[2];
 
  unique_ptr<double[]> sigma2_KJ_block = sigma2_KJ->get_block(sigma2_KJ_id_blocks);
  double* sigma4_ptr = sigma4_IJ_block.get();
  size_t  sigma4_pos_shift = 0 ;
  cout << "[ (sigma4_pos_shifts, Ket_pos_shifts)] = [ " <<  endl;
  //standard sigma from civec algorithm
  for( size_t ii = sigma4_IJ_offsets[0]; ii != iblock_size + sigma4_IJ_offsets[0]; ii++ ) {
    for( size_t jj = sigma4_IJ_offsets[1]; jj != jblock_size + sigma4_IJ_offsets[1]; jj++ ) { 

      double* KKet_ptr = sigma2_KJ_block.get() ;
      size_t  Ket_pos_shift = 0; 
      cout << "[ " ; cout.flush();
      //loops over orbitals in sigma2
      for ( int kk = 0 ; kk!= kblock_size; kk++ ) {                                                                                        
        for ( int ll = 0 ; ll!= lblock_size; ll++ ) { 

          cout << "("  << Ket_pos_shift << ", " <<  sigma4_pos_shift << ") ";
                                                                       
          for ( DetMap iter : KKet_det->phia(ii, jj)) {                                                                                    
                                                                                                                                           
            int KKet_shift = iter.target*lenb - KKet_offset;                                                                               
            if ( (KKet_shift < 0 ) || (KKet_shift >= KKet_block_size) )  { continue; }     
                                                                                           
            int IBra_shift = iter.source*lenb - sigma4_IJ_offsets[4];                  
            if ( (IBra_shift < 0 ) || (IBra_shift >= IBra_block_size) )  { continue; }     
                                                                                                                                           
            size_t loop_limit = lenb > IBra_block_size ? IBra_block_size : lenb;      /*fix */                                             
            loop_limit = KKet_block_size > loop_limit ? loop_limit : KKet_block_size; /*fix */                                             
            double sign = static_cast<double>(iter.sign);                                                                                  
            for( size_t ib = 0; ib != lenb; ib++)                                                                                          
             *(sigma4_ptr+IBra_shift+ib) += (*(KKet_ptr+KKet_shift+ib) * sign);                                                            
          }                                                                                                                                
                                                                                                                                           
          for( size_t ia = 0; ia != lena; ia++) {                                                                                          
            for ( DetMap iter : KKet_det->phib(ii, jj)) {                                                                                  
                                                                                                                                           
              int KKet_shift = iter.target + (ia*lenb) - KKet_offset;                                                                      
              if ( (KKet_shift < 0 ) || (KKet_shift >= KKet_block_size) ) { continue; }    
                                                                                                                                           
              int IBra_shift = iter.source + (ia*lenb) - sigma4_IJ_offsets[4];                                                         
              if ( (IBra_shift < 0 ) || (IBra_shift >= IBra_block_size) ) { continue; }    
                                                                                                                                           
              double sign = static_cast<double>(iter.sign);                                                                                
                                                                                                                                           
              *(sigma4_ptr + IBra_shift) += (*(KKet_ptr+KKet_shift) * sign);                                                               
                                                                                                                                           
            }                                                                                                                              
          } 
          Ket_pos_shift    += KKet_block_size;
          KKet_ptr         += KKet_block_size;                                                                                                     
          sigma4_pos_shift += IBra_block_size;
          sigma4_ptr       += IBra_block_size;                                                                                                    
        }   
      }
      cout <<  " ] " << endl;
    }     
  }
  cout << " got a sigma block .... " ; cout.flush();
  sigma4_IJ->put_block(sigma4_IJ_block, sigma4_IJ_id_blocks ) ; 
  cout << " and put it in the tensor! " << endl;
  return;
 
}       

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// For now just calculate the necessary sigmas and contracts for the 4idx gamma
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_gamma_4idx_tensor(string gamma_4idx_name )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_gamma_4idx_tensor" << endl;

  string sigma_name = "S_"+gamma_4idx_name;

  shared_ptr<GammaInfo> gamma_4idx_info = Gamma_info_map->at(gamma_4idx_name) ;
  shared_ptr<Tensor_<double>> gamma_4idx ;

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

    //contract tensor with vector    

  } else { 
   
    for ( string gamma_2idx_name : gamma_4idx_info->sub_gammas() ) {

      if ( Sigma_data_map->find("S_"+gamma_2idx_name) != Sigma_data_map->end() ) {

        cout << "already got sigma : " "S_"+gamma_2idx_name << endl; 

      } else {
         
        shared_ptr<GammaInfo> gamma_2idx_info =  Gamma_info_map->at(gamma_2idx_name);

        if ( gamma_2idx_info->Bra_info->name() == gamma_2idx_info->Ket_info->name() ) {
           
         build_sigma_2idx_tensor( Gamma_info_map->at(gamma_2idx_name))  ;
 
        } else {
           
           cout << "Must swap Bra and Ket in gamma_info, not implemented yet" << endl;
        }

      }               

    }
 
    shared_ptr<Tensor_<double>> sigma_KijJ = Sigma_data_map->at( "S_"+gamma_4idx_info->sub_gammas(0) );
    shared_ptr<Tensor_<double>> sigma_KklI = Sigma_data_map->at( "S_"+gamma_4idx_info->sub_gammas(1) );

    string IBra_name =  ( Gamma_info_map->at(gamma_4idx_info->sub_gammas(0)) )->Bra_info->name();
    string KBra_name =  ( Gamma_info_map->at(gamma_4idx_info->sub_gammas(1)) )->Bra_info->name();

    shared_ptr<Tensor_<double>> gamma_ij =  Tensor_Calc->contract_tensor_with_vector( sigma_KijJ, CIvec_data_map->at(IBra_name),  make_pair(2,0));
    Gamma_info_map->emplace( "gamma_ij", Gamma_info_map->at(gamma_4idx_info->sub_gammas(0)));
    Gamma_data_map->emplace( "gamma_ij", gamma_ij);
    gamma_2idx_contract_test("gamma_ij");

    shared_ptr<Tensor_<double>> gamma_kl =  Tensor_Calc->contract_tensor_with_vector( sigma_KklI, CIvec_data_map->at(KBra_name),  make_pair(2,0));
    Gamma_info_map->emplace( "gamma_kl", Gamma_info_map->at(gamma_4idx_info->sub_gammas(1)));
    Gamma_data_map->emplace( "gamma_kl", gamma_kl );
    gamma_2idx_contract_test("gamma_kl");
    
    shared_ptr<Tensor_<double>> gamma_4idx_orig_order = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( sigma_KijJ, sigma_KklI, make_pair(2,2) ); 

    shared_ptr<vector<int>> new_order = make_shared<vector<int>>(vector<int> { 1, 0, 2, 3} ); 
    gamma_4idx = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( gamma_4idx_orig_order,  new_order);
    
  }

  cout << " putting " << gamma_4idx_name << " data into map " <<  endl; 
  Gamma_data_map->emplace(gamma_4idx_name , gamma_4idx);

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> 
Gamma_Computer::Gamma_Computer::convert_civec_to_tensor( shared_ptr<const Civec> civector, int state_num ) const {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::convert_civec_to_tensor" << endl;

  //NOTE: must be adapted to handle arbitrary spin sectors
  string civec_name = get_civec_name(state_num, civector->det()->norb(), civector->det()->nelea(), civector->det()->neleb());  
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

#ifndef NDEBUG
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Tests 2idx gamma matrix by taking trace  and checking symmetry. Non-rel only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Gamma_Computer::Gamma_Computer::gamma_2idx_contract_test( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
   bool passed = true;

   cout << endl << "------------------------------------------------------------------------------------------------------" << endl;
   Print_Tensor(Gamma_data_map->at(gamma_name), gamma_name);
   cout << endl << "------------------------------------------------------------------------------------------------------" << endl;
   cout << "Gamma_data_map->at("<<gamma_name<<")->norm() = "<< Gamma_data_map->at(gamma_name)->norm() <<  endl;  
   cout << "Gamma_data_map->at("<<gamma_name<<")->rms()  = "<< Gamma_data_map->at(gamma_name)->rms() <<  endl; 

   shared_ptr<Tensor_<double>> gamma_2idx_trace  =  Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( Gamma_data_map->at(gamma_name),  make_pair(0,1) );
   cout << "------------------------------------------------------------------------------------------------------" << endl; 
   cout << "gamma_2idx_trace->rms() = "<< gamma_2idx_trace->rms() << endl;
   cout << "------------------------------------------------------------------------------------------------------" << endl; 

   shared_ptr<vector<int>>     new_order = make_shared<vector<int>>(vector<int> { 1, 0} ); 

   shared_ptr<Tensor_<double>> gamma_2idx_orig_order = make_shared<Tensor_<double>>( *(Gamma_data_map->at(gamma_name)) );
   shared_ptr<Tensor_<double>> gamma_2idx_transposed = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( gamma_2idx_orig_order, new_order);

   gamma_2idx_transposed->ax_plus_y(-1, gamma_2idx_orig_order) ;


   int nel = Gamma_info_map->at(gamma_name)->Bra_info->nele();
   if ( (abs(gamma_2idx_trace->rms() -nel) > 0.00000001 ) || (gamma_2idx_transposed->rms() > 0.00000001 )  )  passed = false;

   return passed;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Tests 4idx gamma matrix by contracting in different ways and compairing with 2idx gamma matrices.  Only set up for non-rel at present.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Gamma_Computer::Gamma_Computer::gamma_4idx_contract_test( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
   cout << endl << "------------------------------------------------------------------------------------------------------" << endl; 
   Print_Tensor(Gamma_data_map->at(gamma_name));
   cout << endl << "------------------------------------------------------------------------------------------------------" << endl; 
   cout << "Gamma_data_map->at("<<gamma_name<<")->norm() = "<< Gamma_data_map->at(gamma_name)->norm() <<  endl;  
   cout << "Gamma_data_map->at("<<gamma_name<<")->rms()  = "<< Gamma_data_map->at(gamma_name)->rms() <<  endl; 

   shared_ptr<Tensor_<double>> gamma_2idx_from_4idx_A  =  Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( Gamma_data_map->at(gamma_name),  make_pair(2,3) );
   cout << endl << "------------------------------------------------------------------------------------------------------" << endl; 
   Print_Tensor(gamma_2idx_from_4idx_A, "gamma_2idx_from_4idx" ); 
   cout << endl << "------------------------------------------------------------------------------------------------------" << endl; 
   cout << "gamma_2idx_from_4idx_A->norm()  = "<< gamma_2idx_from_4idx_A->norm() << endl;
   cout << "gamma_2idx_from_4idx_A->rms()   = "<< gamma_2idx_from_4idx_A->rms() << endl;
   cout << "------------------------------------------------------------------------------------------------------" << endl; 

   vector<int> ctrs_pos = {0,1,2,3};
   shared_ptr<Tensor_<double>> gamma_2idx_from_4idx_B = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_on_same_tensor( Gamma_data_map->at(gamma_name),  make_pair(0,1) );
   cout << " gamma_2idx_from_4idx_B->norm() = " << gamma_2idx_from_4idx_B->norm()  << endl; 
   cout << " gamma_2idx_from_4idx_B->rms()  = " << gamma_2idx_from_4idx_B->rms()  << endl; 

   gamma_2idx_from_4idx_B->ax_plus_y(-1.0, gamma_2idx_from_4idx_A );
   cout << " gamma_2idx_from_4idx_B->norm() - gamma_2idx_from_4idx_A->norm()  = " << gamma_2idx_from_4idx_B->norm()  << endl; 
   cout << " gamma_2idx_from_4idx_B->rms()  - gamma_2idx_from_4idx_A->rms()   = " << gamma_2idx_from_4idx_B->rms()  << endl; 

   return ( gamma_2idx_from_4idx_B->norm() < 0.00000001 );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Consistency tests for building up gamma 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Gamma_Computer::Gamma_Computer::build_sigma_4idx_tensor_tests(shared_ptr<GammaInfo> gamma_4idx_info ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  bool passed = true;

  shared_ptr<GammaInfo> sigma_2idx_KJ_info = Gamma_info_map->at( gamma_4idx_info->sub_gammas(1) );
  shared_ptr<GammaInfo> sigma_2idx_IK_info = Gamma_info_map->at( gamma_4idx_info->sub_gammas(0) );

  string IJ_Bra_name = gamma_4idx_info->Bra_info->name();
  string IJ_Ket_name = gamma_4idx_info->Ket_info->name();
 
  string IK_Bra_name = sigma_2idx_IK_info->Bra_info->name();
  string IK_Ket_name = sigma_2idx_IK_info->Ket_info->name();

  string KJ_Bra_name = sigma_2idx_KJ_info->Bra_info->name();
  string KJ_Ket_name = sigma_2idx_KJ_info->Ket_info->name();

  if ( IJ_Bra_name != IK_Bra_name ) {
    passed =false;  cout << "IJ_Bra_name != IK_Bra_name   : " <<  IJ_Bra_name << " != " << IK_Bra_name  << endl; 

  } else if ( IJ_Ket_name != KJ_Ket_name ) {
    passed =false;  cout << "IJ_Ket_name != KJ_Ket_name   : " <<  IJ_Ket_name << " != " << KJ_Ket_name  << endl; 

  } else if ( IK_Ket_name != KJ_Bra_name ) {
    passed =false;  cout << "IK_Ket_name != KJ_Bra_name   : " <<  IK_Ket_name << " != " << KJ_Bra_name  << endl; 

  }
 
  return passed;
}


#endif
#endif
// auto in_Bra_range = [&sigma_offsets, &Bra_block_size ]( size_t& pos) {
//     if (( pos > (sigma_offsets[0]+Bra_block_size)) || ( pos < sigma_offsets[0] ) ){
//        cout << "out of sigma range " << endl;
//        cout << "pos            = " << pos << endl;
//        cout << "Bra_offset     = " << sigma_offsets[2] << endl;
//        cout << "Bra_Block_size = " << Bra_block_size << endl;
//        return false;
//     }
//     return true;
//     };
//
//   auto in_Ket_range = [&Ket_offset ]( size_t& pos , size_t Ket_block_size ) {
//      if ( ( pos  > (Ket_offset+ Ket_block_size )) || ( pos < Ket_offset ) ){
//         cout << "out of Ket range " << endl;
//         cout << "pos             = " << pos << endl;
//         cout << "Ket_offset      = " << Ket_offset << endl;
//         cout << "Ket_Block_size  = " << Ket_block_size << endl;
//         return false;
//      }
//      return true;
//      };


