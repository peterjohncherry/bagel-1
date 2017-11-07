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

  cimaxblock = 100; //figure out what is best, maxtile is 10000, so this is chosen to have one index block. Must be consistent if contraction routines are to work...

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
  shared_ptr<vector<string>> gamma_ranges_str  = gamma_info->id_ranges;

  cout << "CIvec_data_map->at("<<Bra_name<<")->norm() = " << CIvec_data_map->at(Bra_name)->norm() << endl; 
  cout << "CIvec_data_map->at("<<Bra_name<<")->rms()  = " << CIvec_data_map->at(Bra_name)->rms() << endl; 

  shared_ptr<vector<bool>> aops = make_shared<vector<bool>>(vector<bool> { true, false } );
  string sigma_name = "S_"+gamma_name;

  if ( Sigma_data_map->find(sigma_name) != Sigma_data_map->end() ){ 
     
     cout << "already got " << sigma_name << ", so use it " << endl;
  
  } else { //should replace this with blockwise and immediate contraction build.

    build_sigma_2idx_tensor( gamma_info );
  
    cout << "sigma tensor = " <<  sigma_name << endl;  
    Print_Tensor(Sigma_data_map->at(sigma_name)); 
    cout << " Sigma_data_map->at("<<sigma_name<<")->rms()  = " <<  Sigma_data_map->at(sigma_name)->rms() << endl;
    cout << " Sigma_data_map->at("<<sigma_name<<")->norm() = " <<  Sigma_data_map->at(sigma_name)->norm() << endl;

    
    cout << "	calculating brabra " << endl;
    double brabra = Tensor_Calc->contract_vectors(  CIvec_data_map->at(Bra_name) , CIvec_data_map->at(Bra_name));
    cout << " brabra  = " << brabra  << endl;

    shared_ptr<Tensor_<double>> gamma_2idx =  Tensor_Calc->contract_tensor_with_vector( Sigma_data_map->at(sigma_name), CIvec_data_map->at(Bra_name),  make_pair(2,0));
    cout << " printing gamma.... " << endl;
    Print_Tensor(gamma_2idx);

 
    Gamma_data_map->emplace( gamma_name, gamma_2idx );

    cout << " Gamma_data_map->at("<< gamma_name<<")->rms()  = " <<  Gamma_data_map->at(gamma_name)->rms() << endl;
    cout << " Gamma_data_map->at("<< gamma_name<<")->norm() = " <<  Gamma_data_map->at(gamma_name)->norm() << endl;

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
 
  // temp hack for aops, must implement SigmaInfo class and map so this is cleaner

  shared_ptr<vector<bool>>   aops = gamma_info->aops; 
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
          
      cout << "block_pos = ["     ; cout.flush() ; for( int  soffset : *block_pos) {cout << soffset << " " ;} cout << "]" << endl;
      cout << "sigma_offsets = [" ; cout.flush() ; for( int  soffset : sigma_offsets) {cout << soffset << " " ;} cout << "]" << endl;
      vector<Index> sigma_id_blocks = *(get_rng_blocks( block_pos, *sigma_ranges));
      build_sigma_block( sigma_tensor, sigma_id_blocks, sigma_offsets, Ket_name ) ;
      cout << "out of buld_sigma_block " << endl; 
    
    } while (fvec_cycle(block_pos, range_lengths, mins ));
    
    Sigma_data_map->emplace(sigma_name, sigma_tensor); 

  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from routines in src/ci/fci/knowles_compute.cc so returns block of a sigma vector in a manner more compatible with
// the Tensor format.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Gamma_Computer::Gamma_Computer::build_sigma_block( shared_ptr<Tensor_<double>> sigma_tens, vector<Index>& sigma_id_blocks,
                                                        vector<int>& sigma_offsets, string Ket_name  ) const {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Gamma_Computer::build_sigma_block" << endl; 

  size_t iblock_size    = sigma_id_blocks[0].size();
  size_t jblock_size    = sigma_id_blocks[1].size();
  size_t Bra_block_size = sigma_id_blocks[2].size();
  size_t Ket_offset = 0;

  const size_t sigma_block_size = Bra_block_size*(iblock_size*norb_+jblock_size);

  cout << " iblock_size      = " << iblock_size << endl;
  cout << " jblock_size      = " << jblock_size << endl;
  cout << " Bra_block_size   = " << Bra_block_size << endl;
  cout << " sigma_block_size = " << sigma_block_size << endl;

  unique_ptr<double[]> sigma_block(new double[sigma_block_size])  ;
  std::fill_n(sigma_block.get(), sigma_block_size, 0.0);

  shared_ptr<const Determinants> Ket_det   = Determinants_map->at(Ket_name);
  shared_ptr<Tensor_<double>>    Ket_Tens  = CIvec_data_map->at(Ket_name);
  shared_ptr<IndexRange>         Ket_range = range_conversion_map->at(Ket_name);
  // size_t lb = Ket_det->lenb() <= Ket_idx_block.size() ? Ket_det_lenb() : Ket_idx_block.size();      
   const int lena = Ket_det->lena();
   const int lenb = Ket_det->lenb();
 
  cout << "iblock_size , jblock_size = " << iblock_size <<" , " << jblock_size << endl;
  // Must be changed to use BLAS, however, you will have to be careful about the ranges; 
  // ci_block needs to have a length which is some integer multiple of lenb.
  for ( Index Ket_idx_block : Ket_range->range()) {  
    cout << " Ket_range = " << Ket_offset  << "-> "  << Ket_offset+Ket_idx_block.size()  << endl;
    cout << " Bra_range = " << sigma_offsets[2]  << "-> "  << sigma_offsets[2] + Bra_block_size << endl;
    vector<Index> Ket_id_vec = {Ket_idx_block};
    unique_ptr<double[]>   Ket_block = make_unique<double[]>(Ket_idx_block.size());
    auto bob   = Ket_Tens->get_block(Ket_id_vec);
     //Ket_block = Ket_Tens->get_block(Ket_id_vec);

    cout << "got ket block " << endl;
    double* Ket_ptr = Ket_block.get();
    size_t Ket_block_size = Ket_idx_block.size();
   
    cout << " sigma_offsets[0],  iblock_size + sigma_offsets[0] = " << sigma_offsets[0] << " , " <<  iblock_size + sigma_offsets[0] << endl;
    cout << " sigma_offsets[1],  jblock_size + sigma_offsets[0] = " << sigma_offsets[1] << " , " <<  jblock_size + sigma_offsets[0] << endl;

    for( size_t ii = sigma_offsets[0]; ii != iblock_size + sigma_offsets[0]; ii++ ) {
      for( size_t jj = sigma_offsets[1]; jj != jblock_size + sigma_offsets[1]; jj++ ) {

        double* sigma_ptr = sigma_block.get() + (ii*iblock_size+jj)*Bra_block_size;

        cout << "sigma alphas ii,jj " << ii << " , " << jj ; cout.flush() ;
        for ( DetMap iter : Ket_det->phia(ii, jj)) {
          if (ii == 4 && jj==4 ) cout << "A1" << endl;
 
          int Ket_shift = iter.target*lenb - Ket_offset;
          if ( (Ket_shift < 0 ) || (Ket_shift >= Ket_block_size) ) {
              cout << "Ket_shift = " << Ket_shift << " out of bounds : Ket_block_size =  " << Ket_block_size   << endl;
              continue;
            } else {
              cout << "Ket_shift = " << Ket_shift << " OK ! " << endl;
            }


          int Bra_shift = iter.source*lenb - sigma_offsets[2];
          if ( (Bra_shift < 0 ) || (Bra_shift >= Bra_block_size) ){
              cout << "Bra_shift = " << Bra_shift << " out of bounds : Bra_block_size =  " << Bra_block_size   << endl;
              continue;
            } else {
              cout << "Bra_shift = " << Bra_shift << " OK ! " << endl;
            }

             
          size_t loop_limit = lenb > Bra_block_size ? Bra_block_size : lenb; //Ket block size and Bra block size should be the same...

          double* sigma_pos = sigma_ptr + Bra_shift;  
          double* Ket_pos = Ket_ptr+Ket_shift;  
          double sign = static_cast<double>(iter.sign);

          if (ii == 4 && jj==4 ) cout << "A1" << endl;
            
          for( size_t ib = 0; ib != loop_limit; ib++, sigma_pos++, Ket_pos++) 
           *sigma_pos += (*Ket_pos * sign);
          
        }

        cout << endl << "     sigma betas " <<  endl;
       
        for( size_t ia = 0; ia != lena; ia++) {
          for ( DetMap iter : Ket_det->phib(ii, jj)) {

            int Ket_shift = iter.target+ia*lenb - Ket_offset;
            if ( (Ket_shift < 0 ) || (Ket_shift >= Ket_block_size) ) {
              cout << "Ket_shift = " << Ket_shift << " out of bounds : Ket_block_size =  " << Ket_block_size   << endl;
              continue;
            } else {
              cout << "Ket_shift = " << Ket_shift << " OK ! " << endl;
            }
            
            int Bra_shift = iter.source+ia*lenb - sigma_offsets[2];
            if ( (Bra_shift < 0 ) || (Bra_shift >= Bra_block_size) ) {
              cout << "Bra_shift = " << Bra_shift << " out of bounds : Bra_block_size =  " << Bra_block_size   << endl;
              continue;
            } else {
              cout << "Bra_shift = " << Bra_shift << " OK ! " << endl;
            }
 
            double sign = static_cast<double>(iter.sign);

            *(sigma_ptr + Bra_shift) += (*(Ket_ptr+Ket_shift) * sign);
            
          }
        } 
        cout << "leaving ii,jj = " << ii << "," << jj << endl;
      }   
    }     
    cout << "leaving orb loop" << endl;
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


