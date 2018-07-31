#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/debugging_utils.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/util/f77.h>
 
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace WickUtils;
using namespace Debugging_Utils;

//#define __DEBUG_TENSOR_ARITHMETIC_UTILS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<size_t> Tensor_Arithmetic_Utils::get_strides_row_major( const std::vector<Index>& block ) { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_strides" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  vector<size_t> stride_vec( block.size() );
  stride_vec.back() = 1 ;
  vector<Index>::const_reverse_iterator b_it = block.crbegin();
  for( vector<size_t>::reverse_iterator sv_it = stride_vec.rbegin()+1; sv_it != stride_vec.rend(); sv_it++ , b_it++ )
    *sv_it = ( *(sv_it -1)) * b_it->size();

  return stride_vec;
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<size_t> Tensor_Arithmetic_Utils::get_strides_column_major( const std::vector<Index>& block ) { 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_strides" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  vector<size_t> stride_vec( block.size() );
  stride_vec.front() = 1 ;
  vector<Index>::const_iterator b_it = block.begin();
  for( vector<size_t>::iterator sv_it = stride_vec.begin()+1; sv_it != stride_vec.end(); sv_it++ , b_it++ )
    *sv_it = ( *(sv_it -1)) * b_it->size();

  return stride_vec;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elems( unique_ptr<double[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems double" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ); };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   double* ptr = some_data.get();
   double total = 0.0;
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += *ptr ;

   return total; 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
complex<double> Tensor_Arithmetic_Utils::sum_elems( unique_ptr<complex<double>[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems complex<double>" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ); };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   complex<double>* ptr = some_data.get();
   complex<double> total = (complex<double>)(0.0);
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += *ptr ;

   return total; 
}
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elem_norms( unique_ptr<double[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems double" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ); };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   return sum_elems( some_data, length ); 
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Tensor_Arithmetic_Utils::sum_elem_norms( unique_ptr<complex<double>[]>& some_data, size_t length  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::sum_elems complex<double>" << endl;
if ( length == 0 ) { throw logic_error ( "this array has no size! Aborting!" ); };  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
   complex<double>* ptr = some_data.get();
   double total = 0.0;
   for ( size_t ii = 0; ii != length ; ii++, ptr++ ) 
     total += abs(*ptr) ;

   return total; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<complex<double>>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<complex<double>>> Tens_in, vector<string>& range_names,
                                                                               shared_ptr<map< string, shared_ptr<IndexRange> > > range_conversion_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 3arg" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> id_ranges(range_names.size());
  for ( int ii = 0; ii != range_names.size() ; ii++ )
    id_ranges[ii] = *range_conversion_map->at(range_names[ii]);

  shared_ptr<Tensor_<complex<double>>> Tens_out = make_shared<Tensor_<complex<double>>>(id_ranges);
  Tens_out->allocate();

  vector<int> range_maxs  = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<complex<double>[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<complex<double>>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<complex<double>>> Tens_in, const vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 2arg" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<complex<double>>> Tens_out = make_shared<Tensor_<complex<double>>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<complex<double>[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<complex<double>>>  Tensor_Arithmetic_Utils::get_sub_tensor_symm( shared_ptr<Tensor_<complex<double>>> Tens_in, const vector<IndexRange>& id_ranges,
                                                                                    const vector<vector<int>>& transforms, const vector<complex<double>>& factors ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor_symm" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<complex<double>>> Tens_out = make_shared<Tensor_<complex<double>>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 
    vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
    if (Tens_in->exists( id_blocks) ) {
       unique_ptr<complex<double>[]> block = Tens_in->get_block( id_blocks );
       Tens_out->put_block( block, id_blocks );

    } else { 
      vector<complex<double>>::const_iterator f_it = factors.begin();
      bool got_symm_partner = false;
      for ( vector<vector<int>>::const_iterator t_it = transforms.begin(); t_it != transforms.end(); t_it++, f_it++ ){
         vector<Index> trans_id_blocks = reorder_vector(*t_it, id_blocks);
         if ( Tens_in->exists( trans_id_blocks ) ) {
           unique_ptr<complex<double>[]> trans_data_block = Tens_in->get_block( trans_id_blocks );
           unique_ptr<complex<double>[]> data_block = Tensor_Arithmetic::Tensor_Arithmetic<complex<double>>::reorder_tensor_data( trans_data_block.get(), *t_it, trans_id_blocks ) ;
           complex<double> factor = *f_it;
           int block_size = Tens_in->get_size(trans_id_blocks);
           int stride_one = 1; 
           zscal_( &block_size, &factor, data_block.get(), &stride_one); 
           Tens_out->put_block( data_block, id_blocks );
           got_symm_partner = true;
           cout << "Obtained tensor block ";  cout.flush(); Debugging_Utils::print_sizes( id_blocks, "" ); cout.flush();
           cout << " from  "; Debugging_Utils::print_sizes( trans_id_blocks,"" ); cout.flush();
           cout << " using symmetry transformation"; print_vector( *t_it ) ; cout << endl;
           break;
         }
         if ( !got_symm_partner ) {
           Debugging_Utils::print_sizes( id_blocks, "id_block with sizes " ); cout.flush(); cout << "not_found" << endl;
           throw logic_error( "Block fetching failed in get sub tensor symm" );
         } 
      }
    }
  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<double>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<double>> Tens_in, vector<string>& range_names,
                                                                      shared_ptr<map< string, shared_ptr<IndexRange> > > range_conversion_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 3arg" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> id_ranges(range_names.size());
  for ( int ii = 0; ii != range_names.size() ; ii++ )
    id_ranges[ii] = *range_conversion_map->at(range_names[ii]);

  shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(id_ranges);
  Tens_out->allocate();

  vector<int> range_maxs = get_range_lengths( id_ranges );
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<double>>  Tensor_Arithmetic_Utils::get_sub_tensor_symm( shared_ptr<Tensor_<double>> Tens_in, const vector<IndexRange>& id_ranges,
                                                                           const vector<vector<int>>& transforms, const vector<double>& factors ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor_symm" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs = get_range_lengths( id_ranges ) ;
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 
    vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
    if ( Tens_in->exists( id_blocks ) ) {
       unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
       Tens_out->put_block( block, id_blocks );

    } else { 
      vector<double>::const_iterator f_it = factors.begin();
      bool got_symm_partner = false;
      for ( vector<vector<int>>::const_iterator t_it = transforms.begin(); t_it != transforms.end(); t_it++, f_it++ ){
         vector<Index> trans_id_blocks = reorder_vector(*t_it, id_blocks );
         if ( Tens_in->exists( trans_id_blocks ) ) {
           unique_ptr<double[]> trans_data_block = Tens_in->get_block( trans_id_blocks );
           unique_ptr<double[]> data_block =  Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_tensor_data( trans_data_block.get(), *t_it, trans_id_blocks ) ;
           int stride_one = 1;
           int block_size = Tens_in->get_size(trans_id_blocks);
           dscal_( block_size , *f_it, data_block.get(), stride_one ); 
           Tens_out->put_block( data_block, id_blocks );
           got_symm_partner = true;
           cout << "Obtained tensor block ";  cout.flush(); Debugging_Utils::print_sizes( id_blocks, "" ); cout.flush();
           cout << " from  "; Debugging_Utils::print_sizes( trans_id_blocks,"" ); cout.flush();
           cout << " using symmetry transformation"; print_vector( *t_it ) ; cout << endl;
           break;
         }
         if ( !got_symm_partner ) {
           Debugging_Utils::print_sizes( id_blocks, "id_block with sizes " ); cout.flush(); cout << "not_found" << endl;
           throw logic_error( "Block fetching failed in get sub tensor symm" );
         } 
      }
    }
  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>
shared_ptr<Tensor_<double>>  Tensor_Arithmetic_Utils::get_sub_tensor( shared_ptr<Tensor_<double>> Tens_in, const vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithemetic_Utils::get_sub_tensor 2arg" << endl; print_sizes( id_ranges, "id_ranges" ); cout << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(id_ranges);
  Tens_out->allocate();
  vector<int> range_maxs  = get_range_lengths( id_ranges );
  vector<int> block_pos(range_maxs.size(),0);  
  vector<int> mins(range_maxs.size(),0);  

  do { 

     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     unique_ptr<double[]> block = Tens_in->get_block( id_blocks );
     Tens_out->put_block( block, id_blocks );

  } while (fvec_cycle_skipper( block_pos, range_maxs, mins ) ); 

  return Tens_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor( shared_ptr<Tensor_<double>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
   cout << "Tensor_Arithmetic_Utils::Print_Tensor " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
  
   if ( Bagel_id_ranges.size() == 1 ) { //lazy switch; should have generic routine 

     assert( Bagel_id_ranges[0].range(0).size() > 0 ); 

     if ( name == "ID" || ( ( Bagel_id_ranges[0].range().size() == 1 ) && ( Bagel_id_ranges[0].range(0).size() == 1 ) ) ) { 
       cout << "Bagel_id_ranges[0].range(0).size() == 1  " << endl; 
       unique_ptr<double[]> data = Tens->get_block( vector<Index> ( 1, Bagel_id_ranges[0].range(0) ) );
       cout << *data.get()  << endl;

     } else { 

       Print_Vector_Tensor_Format( Tens, name );

     }

   } else {

     shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets_sp( Bagel_id_ranges) ;
     
     vector<int> range_lengths  = get_range_lengths( Bagel_id_ranges );
     vector<int> block_pos(range_lengths.size(),0);  
     vector<int> mins(range_lengths.size(),0);  
      
     do {
     
      
       vector<Index> id_blocks = get_rng_blocks( block_pos, Bagel_id_ranges );
      
       if ( !(Tens->exists(id_blocks)) ) { 
     
         cout << "WARNING : TENSOR BLOCK AT " ; cout.flush() ; print_vector( block_pos ) ; cout << " IS NOT STORED FOR PRINTING!! " << endl;
         cout << "NOTE : THIS IS OK IF THERE IS SYMMETRY, BUT OTHERWISE A BIG PROBLEM!! " << endl;
          
       } else {
     
     
         shared_ptr<vector<int>> id_blocks_sizes= make_shared<vector<int>>(id_blocks.size());
         for( int ii = 0 ;  ii != id_blocks.size(); ii++)
           id_blocks_sizes->at(ii) = id_blocks[ii].size();
        
         int id_block_size = Tens->get_size(id_blocks);
       
         unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
         vector<int> Tens_strides = get_Tens_strides_column_major(*id_blocks_sizes);
         double* pos = T_data_block.get();
         int shift = 0 ; 
         if ( block_pos.size() >2 ) { 
           do {
     
             int id_pos_tmp = shift;
             vector<int> id_pos_rel(Tens_strides.size(), 0 );
             for ( int kk = Tens_strides.size()-1 ; kk != 1 ; kk--){
               id_pos_rel[kk] = id_pos_tmp/Tens_strides.at(kk);
               id_pos_tmp -= id_pos_rel[kk]*Tens_strides.at(kk); 
             }

             vector<int> id_pos(block_pos.size());
             for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
               id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii)) + id_pos_rel[ii];

             print_vector( id_pos_rel, "id_pos_rel" ) ; cout << endl;
             print_vector( id_pos, "id_pos" ) ; cout << endl;
     
             for (int jj = 0 ; jj != id_blocks_sizes->at(0) ; jj++ ){
               for (int kk = 0 ; kk != id_blocks_sizes->at(1) ; kk++ ){
                 cout << T_data_block[shift+kk*id_blocks_sizes->at(0)+jj] << " " ; 
               }
               cout <<endl;
             }
     
             shift +=  id_blocks_sizes->at(1)*id_blocks_sizes->at(0);
     
           } while (shift < id_block_size );
     
         } else {
         
           vector<int> id_pos_rel(Tens_strides.size(), 0 );

           vector<int> id_pos(block_pos.size());
           for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
             id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii)) + id_pos_rel[ii];
           
           print_vector( id_pos, "id_pos" ) ; cout << endl;
    
           for (int jj = 0 ; jj != id_blocks_sizes->at(0) ; jj++ ){
             for (int kk = 0 ; kk != id_blocks_sizes->at(1) ; kk++ ){
               cout << T_data_block[shift+kk*id_blocks_sizes->at(0)+jj] << " " ; 
             }
             cout <<endl;
           }
     
           shift += id_blocks_sizes->at(1)*id_blocks_sizes->at(0);
           cout << endl;
         }
       } 
     
     } while ( fvec_cycle_skipper_f2b( block_pos, range_lengths, mins ));
  
   }  

   return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor( shared_ptr<Tensor_<complex<double>>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error( " have not implemented print routines for complex tensors yet!! Aborting!! " ) ;
  return ;
}   

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format( shared_ptr<Tensor_<double>> VecIn, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format" <<endl; 
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  IndexRange VecIn_Idrng = VecIn->indexrange()[0];
  
  for (int ii = 0 ; ii != VecIn_Idrng.range().size(); ii++ ){ 

    vector<Index> id_block = { VecIn_Idrng.range(ii) } ;
    unique_ptr<double[]> VecIn_data_block = VecIn->get_block( id_block );
    double* data_ptr = VecIn_data_block.get();

    for (int jj = 0 ; jj != id_block[0].size() ; jj++ ){ 
      cout << *data_ptr++ << " " ; 
    }   

    cout.flush();
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Vector_Tensor_Format( shared_ptr<Tensor_<complex<double>>> VecIn, string name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error("have not implemented print routines for complex tensors yet!! Aborting !! " ); 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor_row_major( shared_ptr<Tensor_<double>> Tens, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
   cout << "Tensor_Arithmetic_Utils::Print_Tensor " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;

   vector<IndexRange> Bagel_id_ranges = Tens->indexrange();
   shared_ptr<vector<vector<int>>> block_offsets = get_block_offsets_sp( Bagel_id_ranges) ;

   vector<int> range_lengths = get_range_lengths( Bagel_id_ranges ) ;
   vector<int> block_pos (range_lengths.size(),0);  
   vector<int> mins(range_lengths.size(),0);  
    
   do {

     vector<int> id_pos(block_pos.size());
     for ( int ii = 0 ; ii != Bagel_id_ranges.size(); ii++)
       id_pos[ii] = block_offsets->at(ii).at(block_pos.at(ii));
    
     vector<Index> id_blocks = get_rng_blocks( block_pos, Bagel_id_ranges );
    
     if ( Tens->exists(id_blocks) ) { 

       vector<int> id_blocks_sizes(id_blocks.size());

       for( int ii = 0 ;  ii != id_blocks.size(); ii++)
         id_blocks_sizes[ii] = id_blocks[ii].size();
       
       int id_block_size = accumulate(id_blocks_sizes.begin(), id_blocks_sizes.end(), 1 , std::multiplies<int>());


       unique_ptr<double[]>    T_data_block = Tens->get_block(id_blocks);
       vector<int> Tens_strides = get_Tens_strides(id_blocks_sizes);

       for (int ii = 0 ; ii != id_blocks_sizes.size() ; ii++ ) 
         id_blocks_sizes[ii]--;       

       vector<int> rel_id_pos(id_blocks_sizes.size(),0);  
    
       double* pos = T_data_block.get();

       vector<int> rel_id_pos_old(id_blocks_sizes.size()-1, 0);

       do {

          if ( *(rel_id_pos.end()-2) != rel_id_pos_old.back() ){
            rel_id_pos_old.back() = *(rel_id_pos.end()-2);
            cout << endl; 
          }
          
          bool newmat = false;
          for ( int ii = 0; ii != rel_id_pos.size()-2 ; ii++  ){
            if ( rel_id_pos.at(ii) != rel_id_pos_old[ii] ) {
              newmat= true;
            }
          }
          if (newmat ){
            rel_id_pos_old = vector<int>(rel_id_pos.begin(), rel_id_pos.end()-1);
            cout << endl <<  "[ " ; cout.flush();
            for ( int kk = 0 ; kk != rel_id_pos.size() ; kk++ ) { 
              cout << rel_id_pos.at(kk) + block_offsets->at(kk).at(block_pos.at(kk)) << " " ; cout.flush();
            }
            cout << "] " << endl; 
          }   

          cout << *pos++ << " " ; cout.flush();        
 
       } while (fvec_cycle_skipper(rel_id_pos, id_blocks_sizes, mins ) ) ;
 
     } else { 
       
       cout << "WARNING : TENSOR BLOCK AT " ; cout.flush() ; print_vector( block_pos ) ; cout << " IS NOT STORED FOR PRINTING!! " << endl;
       cout << "NOTE : THIS IS OK IF THERE IS SYMMETRY, BUT OTHERWISE A BIG PROBLEM!! " << endl;
        
     }

   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));
 
   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::Print_Tensor_row_major( shared_ptr<Tensor_<complex<double>>> Tens, string name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw runtime_error("have not implemented print routines for complex tensors yet!! Aborting !! " ); 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_Tens_strides( const vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::get_Tens_strides row major" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<int> Tens_strides(range_sizes.size());
  Tens_strides.front() = 1;
  for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
    Tens_strides.at(ii) = Tens_strides.at(ii-1) * range_sizes[ii-1];
   
  return Tens_strides;
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_Tens_strides_column_major( const vector<int>& range_sizes) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::get_Tens_strides_column_major " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<int> Tens_strides(range_sizes.size());
  Tens_strides.front() = 1;
  if (range_sizes.size() > 1 ) {
    for ( int ii = 1  ; ii!= range_sizes.size(); ii++ ) 
      Tens_strides.at(ii) = Tens_strides.at(ii-1) * range_sizes[ii-1];
  }
  return Tens_strides;
 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at back
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::put_ctr_at_back( const vector<int>& orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctr_at_back" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> new_pos = orig_pos;
  if (ctr_pos != orig_pos.back()){
   
   int ii =0;
   for (  ; ii != orig_pos.size(); ii++)
     if ( orig_pos.at(ii) == ctr_pos ) 
       break; 
 
   new_pos.erase(new_pos.begin()+ii);
   new_pos.push_back(ctr_pos); 
 
   if (new_pos.size() != orig_pos.size())
     throw runtime_error("something has gone wrong in index reshuffling");
 
  }
  return new_pos;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//rearranges position vector to have ctr pos at front
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::put_ctr_at_front( const vector<int>& orig_pos , int ctr_pos){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::put_ctr_at_front" <<  endl;  
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (ctr_pos == orig_pos.front()){

    vector<int> new_pos = orig_pos; 
    return new_pos;
 
  } else { 
    vector<int> tmp_pos = orig_pos;
    int ii =0;
    for (  ; ii != orig_pos.size(); ii++)
      if ( tmp_pos[ii] == ctr_pos ) 
        break; 
    
    tmp_pos.erase(tmp_pos.begin()+ii);
    vector<int> new_pos(1, ctr_pos);
    new_pos.insert(new_pos.end(), tmp_pos.begin(), tmp_pos.end() );
   
    assert(new_pos.size() == orig_pos.size());

    return new_pos;
  }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<Index>
Tensor_Arithmetic_Utils::get_rng_blocks( const vector<int>& block_pos, const vector<IndexRange>& id_ranges) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "TensOp_Arithemetic_Utils::get_rng_blocks " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<Index> id_blocks(block_pos.size());
  for( int ii =0; ii != block_pos.size(); ii++ ) 
     id_blocks[ii] = id_ranges[ii].range(block_pos[ii]);
  return id_blocks;
}
////////////////////////////////////////////////////////////////////////////////////////
vector<int>
Tensor_Arithmetic_Utils::get_num_index_blocks_vec(vector<IndexRange>& rngvec) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "TensOp_Arithemetic_Utils::get_num_index_blocks" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> num_id_blocks_vec(rngvec.size()); 
  vector<int>::iterator iter = num_id_blocks_vec.begin();

  for( IndexRange& id_rng : rngvec ){
     *iter = id_rng.range().size()-1;
      ++iter;
  }

  return num_id_blocks_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<size_t>>
Tensor_Arithmetic_Utils::get_sizes_sp(const vector<Index>& Idvec) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_sizes_sp" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  shared_ptr<vector<size_t>> size_vec = make_shared<vector<size_t>>(Idvec.size()); 
  for( int ii = 0 ; ii != Idvec.size() ; ii++ ) 
     (*size_vec)[ii] = Idvec[ii].size();
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
vector<size_t>
Tensor_Arithmetic_Utils::get_sizes(const vector<Index>& index_block) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_sizes" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<size_t> size_vec(index_block.size()); 
  vector<Index>::const_iterator ib_it = index_block.begin();
  for( vector<size_t>::iterator sv_it = size_vec.begin(); sv_it != size_vec.end(); ++sv_it, ++ib_it ) 
     *sv_it = ib_it->size();
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
template<typename T = size_t> 
vector<T>
Tensor_Arithmetic_Utils::get_sizes_m1(const vector<Index>& index_block) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_sizes_m1" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<T> size_vec(index_block.size()); 
  vector<Index>::const_iterator ib_it = index_block.begin();
  for( typename vector<T>::iterator sv_it = size_vec.begin(); sv_it != size_vec.end(); ++sv_it, ++ib_it ) 
     *sv_it = ib_it->size()-1;
  return size_vec;
}
////////////////////////////////////////////////////////////////////////////////////////
template<>
vector<int>
Tensor_Arithmetic_Utils::get_sizes_m1(const vector<Index>& index_block) {
////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_sizes_m1" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<int> size_vec(index_block.size()); 
  vector<Index>::const_iterator ib_it = index_block.begin();
  for( vector<int>::iterator sv_it = size_vec.begin(); sv_it != size_vec.end(); ++sv_it, ++ib_it ) 
     *sv_it = ib_it->size()-1;
  return size_vec;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
size_t
Tensor_Arithmetic_Utils::get_unc_block_size( vector<Index>& idvec, pair<int,int> ctr ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_unc_block_size" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  size_t block_size = 1; 
  for( int ii = 0 ; ii != idvec.size() ; ii++) 
    if ( (ii != ctr.first) && (ii !=ctr.second))
      block_size *= idvec[ii].size();
  
  return  block_size;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
size_t Tensor_Arithmetic_Utils::get_block_size( vector<Index>::iterator startpos, vector<Index>::iterator endpos   ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_block_size" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  size_t block_size = 1; 
  for( vector<Index>::iterator id_it = startpos ; id_it!=endpos; id_it++ ){ 
    block_size *= id_it->size();
  }
  return  block_size;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int> Tensor_Arithmetic_Utils::get_range_lengths(const vector<IndexRange>& indexranges ) { 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS_VERBOSE
cout << "Tensor_Arithmetic_Utils::get_range_lengths" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> range_lengths(indexranges.size() ); 
  vector<IndexRange>::const_iterator ir_it = indexranges.begin();
  for (vector<int>::iterator rl_it = range_lengths.begin(); rl_it!= range_lengths.end(); rl_it++, ir_it++ )
    *rl_it = ir_it->range().size()-1; 

  return range_lengths;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Checks all the index ranges in idx_pos have the same length and number of index blocks
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::check_contracted_indexes( std::vector<IndexRange>&  idx_block, std::vector<int>& idx_pos ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::check_contracted_indexes " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  size_t block_len = idx_block[idx_pos.front()].size();
  for ( int ii = 1; ii != idx_pos.size() ; ii++ ) 
    if (block_len != idx_block[idx_pos[ii]].size()){
      cout << "ranges of contracted indexes " << idx_pos[ii] << " and " <<  idx_pos[ii-1]  << "have different lengths" << endl;
      assert(block_len == idx_block[idx_pos[ii]].size());
    }

  block_len = idx_block[idx_pos.front()].range().size();
  for ( int ii = 1; ii != idx_pos.size() ; ii++ ) 
    if (block_len != idx_block[idx_pos[ii]].range().size()){
      cout << "ranges of contracted indexes " << idx_pos[ii] << " and " <<  idx_pos[ii-1]  << "have same lengths but different number of index blocks";
      assert(block_len == idx_block[idx_pos[ii]].range().size());
    }

  return;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////
// Outputs the interval [start_orb, end_orb) of the index blocks at block_pos
// Note the interval is relative to the input IndexRanges.
////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<size_t, size_t>>>
Tensor_Arithmetic_Utils::get_block_start( shared_ptr<vector<IndexRange>> id_ranges, 
                                          shared_ptr<vector<int>> block_pos         ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "TensOp_Arithmetic_Utils::get_block_start" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<pair<size_t,size_t>> block_start_end(block_pos->size());
  for (int ii = 0 ; ii != block_start_end.size() ; ii++){
    size_t block_start = 0;

    for (int jj = 0 ; jj != block_pos->at(ii); jj++) 
      block_start += id_ranges->at(ii).range(jj).size();

    block_start_end[ii] = make_pair(block_start, block_start+id_ranges->at(ii).range(block_pos->at(ii)).size());

  }

  return make_shared<vector<pair<size_t,size_t>>>(block_start_end);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<vector<vector<size_t>>,vector<vector<size_t>>> Tensor_Arithmetic_Utils::get_block_start_ends( const vector<IndexRange>& ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
 cout << "Tensor_Arithmetic_Utils::get_block_start_ends" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<vector<size_t>> block_starts(ranges.size());
  vector<vector<size_t>> block_ends(ranges.size());
   
  for ( int  ii = 0 ; ii != ranges.size() ; ii++ ) { 

    vector<Index> blocks = ranges[ii].range();
    vector<size_t> cml_sizes(blocks.size()); 
    cml_sizes[0] = blocks.front().size()-1;

    for ( int  jj = 1 ; jj != blocks.size(); jj++ ) 
      cml_sizes[jj] = cml_sizes[jj-1]+blocks[jj].size(); 

 
    vector<size_t> tmp_starts(cml_sizes.size());
    tmp_starts.front() = 0;

    for ( int qq = 1 ; qq != cml_sizes.size() ; qq++ ) 
      tmp_starts[qq] = cml_sizes[qq-1]+1;          

    block_ends[ii] =  cml_sizes;
    block_starts[ii] = tmp_starts;

  }

  for ( int qq = 0 ; qq != block_starts.size() ; qq++ ) { 
    cout << "block_starts["<< qq <<"] = "; cout.flush();  print_vector( block_starts[qq], ""); cout << endl;
    cout << "block_ends["<< qq <<"] = "; cout.flush();  print_vector( block_ends[qq], ""); cout << endl;
  }

  return make_pair( block_starts, block_ends ); 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<vector<int>>> Tensor_Arithmetic_Utils::get_block_offsets_sp( const vector<IndexRange>&  ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
 cout << "Tensor_Arithmetic_Utils::get_block_offsets_sp" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<vector<vector<int>>> block_offsets = make_shared<vector<vector<int>>>(ranges.size());
   
  for ( int  ii = 0 ; ii != ranges.size() ; ii++ ) { 
    vector<int> block_offset(ranges[ii].range().size()); 
    block_offset[0] = 0;
    for ( int  jj = 1 ; jj != ranges[ii].range().size() ; jj++ ) { 
      block_offset[jj] = block_offset[jj-1]+ranges[ii].range(jj).size(); 
    }    
    block_offsets->at(ii) =  block_offset;
  }
  
  return block_offsets ; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<vector<int>> Tensor_Arithmetic_Utils::get_block_offsets( const vector<IndexRange>&  ranges ) { 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
 cout << "Tensor_Arithmetic_Utils::get_block_offsets" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  vector<vector<int>> block_offsets(ranges.size());
  vector<vector<int>>::iterator bos_it = block_offsets.begin();
   
  for ( vector<IndexRange>::const_iterator r_it = ranges.begin(); r_it != ranges.end(); ++r_it, ++bos_it ){
    vector<int> block_offset(r_it->range().size()); 
    block_offset.front() = 0;
    for ( int  jj = 1 ; jj != r_it->range().size() ; jj++ ) { 
      block_offset[jj] = block_offset[jj-1]+r_it->range(jj).size(); 
    }    
    *bos_it =  block_offset;
  }
  
  return block_offsets ; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front, reverse vector version. Note, this will change the input vector orig_pos. Do not change ctrs_pos. 
// The reversal so that when the indices are flattened matrix multiplication can be used.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_reversed_ctrs_at_front( vector<int>& ids_pos, vector<int> ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_reversed_ctrs_at_front" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }

  reverse(ctrs_pos.begin(), ctrs_pos.end());
  ids_pos.insert(ids_pos.begin(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_front( vector<int>& ids_pos, vector<int>& ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_front" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }
  ids_pos.insert(ids_pos.begin(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at front,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_front( vector<int>& ids_pos, pair<int,int>& ctr_pair) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_front pair" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   vector<int> ctrs_pos = {ctr_pair.first, ctr_pair.second};

   put_ctrs_at_front(ids_pos, ctrs_pos );
 
   return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put ctrs at back,  vector version. Note, this will change the input vector orig_pos
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Tensor_Arithmetic_Utils::put_ctrs_at_back( vector<int>& ids_pos, vector<int>& ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::put_ctrs_at_back" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for ( int ctr_pos : ctrs_pos) 
    for ( vector<int>::iterator iter = ids_pos.begin(); iter != ids_pos.end(); iter++  ) 
      if (  *iter == ctr_pos  ){
        ids_pos.erase(iter); 
        break;
      }

  ids_pos.insert(ids_pos.end(), ctrs_pos.begin(), ctrs_pos.end());

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<IndexRange> Tensor_Arithmetic_Utils::get_indexrange_from_dimension_vector( vector<size_t> dimensions, vector<size_t> maxblock  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
  cout << "Tensor_Arithmetic_Utils::get_indexrange_from_dimension_vector" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> index_ranges(dimensions.size());
  for ( int ii = 0 ; ii != dimensions.size() ; ii++ ) 
    index_ranges[ii] = IndexRange( dimensions[ii], maxblock[ii]);

  return index_ranges;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::print_tensor_with_indexes( shared_ptr<Tensor_<double>> Tens, string name , bool print_zero_blocks ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::print_tensor_with_indexes " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;
   cout.precision(12);

   vector<IndexRange> id_ranges = Tens->indexrange();
  
   if ( id_ranges.size() == 1 ) { 

     assert( id_ranges[0].range(0).size() > 0 ); 

     if ( name == "ID" || ( ( id_ranges[0].range().size() == 1 ) && ( id_ranges[0].range(0).size() == 1 ) ) ) { 
       cout << "id_ranges[0].range(0).size() == 1  " << endl; 
       unique_ptr<double[]> data = Tens->get_block( vector<Index> ( 1, id_ranges[0].range(0) ) );
       cout << *data.get()  << endl;

     } else {
       size_t elem_num = 0;
       auto compare_abs = []( int aa, int bb ) { return (abs(aa) < abs(bb)); };
       for ( const Index& block : id_ranges[0].range() ){
         unique_ptr<double[]> data = Tens->get_block( vector<Index> ( 1, block  ) );
         if ( print_zero_blocks || ( abs(*(max_element(data.get(), data.get()+block.size(), compare_abs )))  > ( 1.0e-14 ) ) ) { 
           size_t endpos = elem_num+block.size();
           for ( double* d_it = data.get() ; elem_num !=  endpos; ++elem_num, ++d_it )
             cout << "[ " << elem_num << " ] : " << *d_it << endl;
         } else {
           cout << endl <<  "-------- elems from "<<  elem_num << " to "; cout.flush();
           elem_num += block.size();
           cout << elem_num << " are zero; not printing ----------------" << endl << endl;
         }
       }
     } cout <<endl;

   } else {

     int t_order = id_ranges.size();

     vector<vector<size_t>> block_starts;
     vector<vector<size_t>> block_ends;
     {
       pair<vector<vector<size_t>>,vector<vector<size_t>>> block_start_ends =  get_block_start_ends( id_ranges );
       block_starts = block_start_ends.first;
       block_ends = block_start_ends.second;
     }

     vector<int> elem_pos_maxs( t_order );
     vector<int> elem_pos_mins( t_order,0);  
     vector<int> elem_pos( t_order,0);  
      
     vector<int> block_pos_maxs = get_range_lengths( id_ranges );
     vector<int> block_pos_mins( t_order, 0);  
     vector<int> block_pos(t_order, 0);  
    
     //TODO will not work if set to bool, find out why 
     auto compare_abs = []( int aa, int bb ) { return  (abs(aa) < abs(bb)) ; };
 
     do {
      
       vector<int> block_mins(t_order,0);     
       vector<int> rel_elem_pos(t_order,0);     
       
       vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
       assert(Tens->exists(id_blocks));
      
       vector<int> block_maxs(t_order);
       { 
       vector<Index>::iterator ib_it = id_blocks.begin();
       for ( vector<int>::iterator bm_it = block_maxs.begin(); bm_it != block_maxs.end() ; bm_it++, ib_it++ )
         *bm_it = (ib_it->size())-1;
       }

       {
       vector<int>::iterator epmn_it = elem_pos_mins.begin();
       vector<int>::iterator epmx_it = elem_pos_maxs.begin();
       vector<vector<size_t>>::iterator bs_it =  block_starts.begin();
       vector<vector<size_t>>::iterator be_it =  block_ends.begin();
       for ( vector<int>::iterator bp_it = block_pos.begin(); bp_it != block_pos.end(); bp_it++, epmn_it++, epmx_it++, bs_it++, be_it++ ) {
         *epmn_it = (*bs_it)[*bp_it];
         *epmx_it = (*be_it)[*bp_it];
       }
       }
       elem_pos = elem_pos_mins;

       print_vector(elem_pos_mins, "elem_pos_mins" ); cout.flush(); print_vector(elem_pos_maxs, "   elem_pos_maxs" ); cout.flush();
       print_vector(block_pos, "block_pos"); cout.flush(); print_sizes(id_blocks, "   id_blocks" );
       unique_ptr<double[]> data_block = Tens->get_block( id_blocks );
       double* ptr = data_block.get();
       if ( print_zero_blocks || ( abs(*(max_element(data_block.get(), data_block.get()+Tens->get_size(id_blocks), compare_abs ))) > ( 1.0e-14 ) ) ) { 
         cout << endl;
         do {  
           print_vector( elem_pos , "" ); cout << " : ";cout.flush(); cout << *ptr << endl;
           ++ptr;
           fvec_cycle_skipper_f2b( elem_pos, elem_pos_maxs, elem_pos_mins);
         } while( fvec_cycle_skipper_f2b(rel_elem_pos, block_maxs, block_mins) );
         ptr -= Tens->get_size( id_blocks );

       } else { 
         cout << " !!! No non-zero elements; not printing !!! " << endl;
         do {  
           fvec_cycle_skipper_f2b( elem_pos, elem_pos_maxs, elem_pos_mins);
         } while( fvec_cycle_skipper_f2b(rel_elem_pos, block_maxs, block_mins) );
       } 
       
     } while ( fvec_cycle_skipper_f2b( block_pos, block_pos_maxs, block_pos_mins ));
   }  

   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::print_tensor_with_indexes( shared_ptr<Tensor_<std::complex<double>>> Tens, string name , bool print_zero_blocks ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::print_tensor_with_indexes " << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "---------------------------- " << name <<  " ----------------------------" << endl;
   throw logic_error( " no printing for complex tensors yet" ); 
   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::set_test_elems( shared_ptr<Tensor_<double>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::set_test_elems" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   vector<IndexRange> id_ranges = Tens->indexrange();
  
   vector<vector<size_t>> block_starts;
   vector<vector<size_t>> block_ends;
   {
     pair<vector<vector<size_t>>,vector<vector<size_t>>> block_start_ends =  get_block_start_ends( id_ranges );
     block_starts = block_start_ends.first;
     block_ends = block_start_ends.second;
   }

   vector<double> factor_vec( id_ranges.size());
   double exponent = 0.0;
   for ( vector<double>::reverse_iterator fv_it = factor_vec.rbegin(); fv_it != factor_vec.rend(); fv_it++, exponent+=1.0 ) 
     *fv_it = pow( 10.0, exponent ); 
   
   print_vector( factor_vec, "factor_vec"); cout << endl;

   auto dot_vector = [&factor_vec]( vector<int>& elem_pos_vec ) { 
                     double elem = 0.0;
                     vector<int>::iterator epv_it = elem_pos_vec.begin(); 
                     for ( vector<double>::iterator fv_it = factor_vec.begin(); fv_it != factor_vec.end(); fv_it++, epv_it++ )
                       elem += (*fv_it) * (*epv_it);
                     return elem;
                     };
 
   int t_order = id_ranges.size();

   vector<int> elem_pos_maxs( t_order );
   vector<int> elem_pos_mins( t_order );  
   vector<int> elem_pos( t_order,0);  
    
   vector<int> block_pos_maxs = get_range_lengths( id_ranges );
   vector<int> block_pos_mins( t_order, 0);  
   vector<int> block_pos(t_order, 0);  
    
   do {
    
     vector<int> block_mins(t_order,0);     
     vector<int> rel_elem_pos(t_order,0);     
     
     vector<Index> id_blocks = get_rng_blocks( block_pos, id_ranges );
     assert(Tens->exists(id_blocks));
    
     vector<int> block_maxs(t_order);
     { 
     vector<Index>::iterator ib_it = id_blocks.begin();
     for ( vector<int>::iterator bm_it = block_maxs.begin(); bm_it != block_maxs.end() ; bm_it++, ib_it++ )
       *bm_it = (ib_it->size())-1;
     }

     {
     vector<int>::iterator epmn_it = elem_pos_mins.begin();
     vector<int>::iterator epmx_it = elem_pos_maxs.begin();
     vector<vector<size_t>>::iterator bs_it =  block_starts.begin();
     vector<vector<size_t>>::iterator be_it =  block_ends.begin();
     for ( vector<int>::iterator bp_it = block_pos.begin(); bp_it != block_pos.end(); bp_it++, epmn_it++, epmx_it++, bs_it++, be_it++ ) {
       *epmn_it = (*bs_it)[*bp_it];
       *epmx_it = (*be_it)[*bp_it];
     }
     }
     elem_pos = elem_pos_mins;

     unique_ptr<double[]> data_block = Tens->get_block( id_blocks );
     double* ptr = data_block.get();
     print_vector(block_pos, "block_pos"); cout << endl;
     do {  
       *ptr = dot_vector( elem_pos );
       ++ptr;
       fvec_cycle_skipper_f2b( elem_pos, elem_pos_maxs, elem_pos_mins );
     } while( fvec_cycle_skipper_f2b(rel_elem_pos, block_maxs, block_mins) );
     Tens->put_block( data_block, id_blocks);
   } while ( fvec_cycle_skipper_f2b( block_pos, block_pos_maxs, block_pos_mins ));

   return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<>     
void Tensor_Arithmetic_Utils::set_test_elems( shared_ptr<Tensor_<complex<double>>> Tens, string name  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOR_ARITHMETIC_UTILS
cout << "Tensor_Arithmetic_Utils::get_test_tensor" << endl;
#endif /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  throw logic_error ( "cannot set complex test tensor" ); 
  return;
}  

