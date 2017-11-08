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
Gamma_Computer::Gamma_Computer::Gamma_Computer( shared_ptr< map< string, shared_ptr<GammaInfo>>>          GammaMap_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    CIvec_data_map_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    Sigma_data_map_in,
                                                shared_ptr< map< string, shared_ptr<Tensor_<double>>>>    Gamma_data_map_in,
                                                shared_ptr< map< string, shared_ptr<const Determinants>>> Determinants_map_in,
						shared_ptr< map< string, shared_ptr<IndexRange>>>         range_conversion_map_in ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::Gamma_Computer" << endl;
  maxtile  = 10000;

  GammaMap             = GammaMap_in;             cout << "set GammaMap            "<< endl; 
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
void Gamma_Computer::Gamma_Computer::get_gamma_tensor( string gamma_name ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::Gamma_Computer::get_gammas"  << endl;
  cout << "Gamma name = " << gamma_name << endl;

  if( Gamma_data_map->find(gamma_name) != Gamma_data_map->end()){ 

    cout << "already have data for " << gamma_name << endl;  
 
  } else { 
    
    //for now just use specialized routines, this must be made generic at some point
    if (GammaMap->at(gamma_name)->id_ranges->size() == 2 ) { 
      build_gamma_2idx_tensor( gamma_name ) ;
      cout << "------------------ "<<  gamma_name  << " ---------------------" << endl; 
      Print_Tensor(Gamma_data_map->at(gamma_name));
      cout << "Gamma_data_map->at("<<gamma_name<<")->norm() = "<< Gamma_data_map->at(gamma_name)->norm() <<  endl; 
      cout << "Gamma_data_map->at("<<gamma_name<<")->rms()  = "<< Gamma_data_map->at(gamma_name)->rms() <<  endl; 
      
    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 4 ) { 
       cout << " 4-index stuff not implemented, cannot calculate " << gamma_name << endl;

    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 6 ) { 
       cout << " 6-index stuff not implemented, cannot calculate " << gamma_name << endl;

    } else if (GammaMap->at(gamma_name)->id_ranges->size() == 8 ) { 
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

  shared_ptr<GammaInfo> gamma_info      =  GammaMap->at(gamma_name);

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
 
  shared_ptr<vector<string>> orb_ranges_str = gamma_info->id_ranges; 

  string sigma_name = "S_"+gamma_info->name;

  string Bra_name = gamma_info->Bra_info->name();
  string Ket_name = gamma_info->Ket_info->name();

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

  } else { 
    
    shared_ptr<Tensor_<double>>  Bra_civec = CIvec_data_map->at(Bra_name);
    shared_ptr<const Determinants> Bra_det = Determinants_map->at(Bra_name);
    
    shared_ptr<vector<IndexRange>> orb_ranges  = Get_Bagel_IndexRanges( orb_ranges_str );   
    
    auto sigma_ranges = make_shared<vector<IndexRange>>(vector<IndexRange> { orb_ranges->at(0), orb_ranges->at(1), Bra_civec->indexrange()[0]} );

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
// For now just calculate the necessary sigmas and contracts for the 4idx gamma
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
Gamma_Computer::Gamma_Computer::build_gamma_4idx_tensor(shared_ptr<GammaInfo> gamma_4idx_info )  {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "build_sigma_4idx_tensor" << endl;

  string sigma_name = "S_"+gamma_4idx_info->name;

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
   
    cout << "already got " << sigma_name << endl;

  } else { 
   
  for ( string gamma2_name : gamma_4idx_info->sub_gammas() ) {

        if ( Sigma_data_map->find("S_"+gamma2_name) != Sigma_data_map->end() ) {

        cout << "already got sigma : " "S_"+gamma2_name << endl; 

      } else {
         
        shared_ptr<GammaInfo> gamma_2idx_info =  GammaMap->at(gamma2_name);

        if ( gamma_2idx_info->Bra_info->name() == gamma_2idx_info->Ket_info->name() ) {
     
         build_sigma_2idx_tensor( GammaMap->at(gamma2_name))  ;

         shared_ptr<Tensor_<double>> sigma_KijJ = Sigma_data_map->at( "S_"+gamma_4idx_info->sub_gammas(0) );
         shared_ptr<Tensor_<double>> sigma_KklI = Sigma_data_map->at( "S_"+gamma_4idx_info->sub_gammas(1) );

         shared_ptr<vector<int>> new_order = make_shared<vector<int>>(vector<int> { 1, 0, 2, 3} ); 
         
         shared_ptr<Tensor_<double>> gamma_4idx_orig_order = Tensor_Arithmetic::Tensor_Arithmetic<double>::contract_different_tensors( sigma_KijJ, sigma_KklI, make_pair(2,2) ); 
         shared_ptr<Tensor_<double>> gamma_4idx = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( gamma_4idx_orig_order,  new_order);
         
         Gamma_data_map->emplace(gamma_4idx_info->name , gamma_4idx);

       
        } else {
           
           cout << "Must swap Bra and Ket in gamma_info, not implemented yet" << endl;
        }
 
      }               

    }



  }

  return;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
// TODO Skipping for CI blocks breaks in certain cases; orb block serperation is OK. 
// Left in for now as innocuous  unless CI vector becomes very large, and outline is helpful.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigma_block( shared_ptr<Tensor_<double>> sigma_tens, vector<Index>& sigma_id_blocks,
                                                        vector<int>& sigma_offsets, string Ket_name  ) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigma_block" << endl; 

  size_t iblock_size    = sigma_id_blocks[0].size();
  size_t jblock_size    = sigma_id_blocks[1].size();
  size_t Bra_block_size = sigma_id_blocks[2].size();
  size_t Ket_offset = 0;

  const size_t sigma_block_size = Bra_block_size*iblock_size*jblock_size;

  unique_ptr<double[]> sigma_block(new double[sigma_block_size])  ;
  std::fill_n(sigma_block.get(), sigma_block_size, 0.0);

  shared_ptr<const Determinants> Ket_det   = Determinants_map->at(Ket_name);
  shared_ptr<Tensor_<double>>    Ket_Tens  = CIvec_data_map->at(Ket_name);
  shared_ptr<IndexRange>         Ket_range = range_conversion_map->at(Ket_name);

  const int lena = Ket_det->lena();
  const int lenb = Ket_det->lenb();
  double* sigma_block_test;
  
  sigma_block_test = (double*) malloc(sizeof(double)*(sigma_block_size)); 

  // Must be changed to use BLAS, however, you will have to be careful about the ranges; 
  // ci_block needs to have a length which is some integer multiple of lenb.
  for ( Index Ket_idx_block : Ket_range->range()) {  

    unique_ptr<double[]>  Ket_block = Ket_Tens->get_block( vector<Index> {Ket_idx_block} );

    double* Ket_ptr = Ket_block.get();
    size_t Ket_block_size = Ket_idx_block.size();

    for( size_t ii = sigma_offsets[0]; ii != iblock_size + sigma_offsets[0]; ii++ ) {
      for( size_t jj = sigma_offsets[1]; jj != jblock_size + sigma_offsets[1]; jj++ ) {

        double* sigma_ptr = sigma_block.get() + (ii*iblock_size+jj)*Bra_block_size;
        
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
      }   
    }     
    Ket_offset += Ket_idx_block.size(); 
  }
  cout << " got a sigma block .... " ; cout.flush();
  sigma_tens->put_block(sigma_block, sigma_id_blocks ) ; 
  cout << " and put it in the tensor! " << endl;
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


