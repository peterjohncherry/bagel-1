#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_orb_exc_deriv.h>

using namespace std;
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::add_gamma( const shared_ptr<Range_Block_Info> block_info, shared_ptr<vector<bool>> trans_aops  ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "void GammaGenerator_OrbExcDeriv::add_gamma " << endl;

  block_idxs_ = vector<string>( std_ids_->size() );
  idxs_trans_ = block_info->idxs_trans();
  block_aops_ = trans_aops;
  block_aops_rngs_ = block_info->orig_rngs_ch(); // May be transformed by Bra or Ket.

  block_rngs_ = make_shared<vector<string>>( block_info->orig_rngs()->size() ) ; // After the below definition, this remains unchanged.
  vector<string>::iterator br_it = block_rngs_->begin();
  for ( vector<int>::iterator it_it = idxs_trans_->begin() ; it_it != idxs_trans_->end(); it_it++ , br_it++ )
    *br_it = (*(block_info->unique_block_->orig_rngs_))[*it_it];

  shared_ptr<vector<int>>  idxs_trans_inverse_ = block_info->idxs_trans_inverse();
  std_rngs_ = *(block_info->unique_block_->orig_rngs_);
  standard_order_  = *(block_info->idxs_trans());

  cout << "block_info->name = " << block_info->name() << endl;
  cout << "block_info->full_op_name = " << block_info->full_op_name() << endl;
  cout << "target_op_= " << target_op_ << endl;

  // TODO change; this will be problematic when we have multiple operators with the same name but different states.
  target_block_start_ = 0 ;
  target_block_end_ = 0 ;
  for ( vector<shared_ptr<Range_Block_Info>>::iterator rb_it = block_info->range_blocks()->begin() ; rb_it != block_info->range_blocks()->end(); rb_it++ ) {
    target_block_end_ += (*rb_it)->num_idxs_; 
    if ( (*rb_it)->full_op_name()[0] == target_op_[0] )
      target_block_name_ =  (*rb_it)->full_op_name(); 
      
    target_block_start_ += (*rb_it)->num_idxs_; 
  }

  cout << "target_block_name_ = " ; cout.flush(); cout << target_block_name_  << endl;
 
  {
  vector<int>::iterator it_it = idxs_trans_->begin();
  for ( vector<string>::iterator bi_it = block_idxs_.begin(); bi_it != block_idxs_.end(); bi_it++, it_it++ )
    *bi_it = (*std_ids_)[ *it_it ];

  }
  print_vector( block_idxs_, "block_idxs" ); cout <<endl;

  cout << endl;
  cout << "----------- gamma orb deriv def -------------------" << endl;
  print_vector( std_rngs_ ,    " unique_block_ "); cout <<endl;
  print_vector(*idxs_trans_ ,  " idxs_trans_   "); cout << endl;
  print_vector(*block_rngs_ ,  " block_rngs    "); cout <<endl;
  cout << endl;

  int ii = 0 ;
  block_to_std_order_ = vector<int>(idxs_trans_->size());
  for ( vector<int>::iterator so_it = idxs_trans_->begin() ; so_it != idxs_trans_->end() ; ++so_it, ++ii ) 
    block_to_std_order_[*so_it] = (ii);

  shared_ptr<vector<int>> ids_pos = make_shared<vector<int>>( std_rngs_.size() );
  iota( ids_pos->begin(), ids_pos->end(), 0 );

  shared_ptr<vector<pair<int,int>>> A_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> target_target_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> target_A_deltas_pos = make_shared<vector<pair<int,int>>>(0);
  
  pair< double, double >  factors = block_info->factors();
  cout << "factors = (" << factors.first << "," << factors.second << ")" << endl;

  //TODO  Change to specialized class
  gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediate_Base >>>();
  gamma_vec_->push_back(make_shared<GammaIntermediate_OrbExcDeriv<DataType>>( ids_pos, A_deltas_pos, target_A_deltas_pos, target_target_deltas_pos, factors ));
  final_gamma_vec_ = make_shared<vector<shared_ptr<GammaIntermediate_Base >>>(0);
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cannot pass shared_ptr to gammaintermediate, as push back can potentially result in the vector being moved,
// which messes up the pointer inside the shared_ptr.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::swap( int ii, int jj, int kk ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "GammaGenerator_OrbExcDeriv<DataType>::swap ii = " << ii << " jj = " << jj << " kk = " << kk << endl;

  shared_ptr<GammaIntermediate_Base> gint =  gamma_vec_->at(kk);

  int j_pos = (*((gint)->ids_pos_))[jj];
  int i_pos = (*((gint)->ids_pos_))[ii];

  (*((gint)->ids_pos_))[ii] = j_pos;
  (*((gint)->ids_pos_))[jj] = i_pos;

  if ( ( (*block_aops_rngs_)[ j_pos ] == (*block_aops_rngs_)[ i_pos ]) &&( (*block_aops_)[ i_pos ] != (*block_aops_)[ j_pos ]) ){

      shared_ptr<vector<int>> new_ids_pos = make_shared<vector<int>>( (gint)->ids_pos_->size()-2);
      vector<int>::iterator nip_it = new_ids_pos->begin();
      vector<int>::iterator gip_it = gint->ids_pos_->begin();
      for( int qq = 0 ; qq != gint->ids_pos_->size() ; qq++, gip_it++) {
        if ( (qq != ii) && (qq != jj)){
         *nip_it = *gip_it;
          ++nip_it;
        }
      }

      shared_ptr<pint_vec> new_target_target_deltas_pos = make_shared<pint_vec>( *(gint->target_target_deltas_pos()) );
      shared_ptr<pint_vec> new_target_A_deltas_pos = make_shared<pint_vec>( *(gint->target_A_deltas_pos()) );
      shared_ptr<pint_vec> new_A_A_deltas_pos = make_shared<pint_vec>( *(gint->deltas_pos_) );

      pair<int,int> new_delta = (*block_aops_)[ j_pos ] ? make_pair( j_pos, i_pos ) : make_pair( i_pos , j_pos);

      
      if ( (j_pos >= target_block_start_) && ( j_pos < target_block_end_ )  ) {

        if ( (i_pos >= target_block_start_) &&  ( i_pos < target_block_end_ )  ) {
          cout << "T x T" << endl;
	  new_target_target_deltas_pos->push_back(new_delta);
        } else {
          cout << "T x A" << endl;
	  new_target_A_deltas_pos->push_back(new_delta);
        }
      } else if ( i_pos >= target_block_start_ && i_pos < target_block_end_ ) { 
        cout << "T x A" << endl;
        new_target_A_deltas_pos->push_back(new_delta);
      } else {  
        cout << "A x A" << endl;
        new_A_A_deltas_pos->push_back(new_delta);
      } 
  
    shared_ptr<GammaIntermediate_OrbExcDeriv<DataType>> new_gamma = 
         make_shared<GammaIntermediate_OrbExcDeriv<DataType>>( new_ids_pos, new_A_A_deltas_pos, new_target_A_deltas_pos, new_target_target_deltas_pos, (gint)->factors_ );

    gamma_vec_->push_back(new_gamma);

  }

  //Note that these are factors for the real and imaginary part, they are _not_ the real and imaginary part of the factor
  pair<double,double> anti_herm_fac = make_pair( gint->factors_.first*-1.0 ,  gint->factors_.second*-1 );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::add_Acontrib_to_map( int kk, string bra_name, string ket_name ){  // e.g. ++++----
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "GammaGenerator_OrbExcDeriv<DataType>::add_Acontrib_to_map" << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void GammaGenerator_OrbExcDeriv<DataType>::transformation_tester( int kk  ){  // e.g. ++++----
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "GammaGenerator_OrbExcDeriv<DataType>::add_Acontrib_to_map" << endl;


  shared_ptr<GammaIntermediate_Base> gint = gamma_vec_->at(kk);
  shared_ptr<vector<pair<int,int>>> deltas_pos = gint->deltas_pos_;

  cout << endl << endl << endl;
  cout << "----------- gamma orb deriv Contribs -------------------" << endl;
  print_vector( *(gint->ids_pos_) ,                         " ids_pos     _            "); cout <<endl;
  print_pair_vector( *( gint->deltas_pos_),                 " A_A_deltas_pos           "); cout << endl;
  print_pair_vector( *( gint->target_target_deltas_pos()),  " target_target_deltas_pos "); cout << endl;
  print_pair_vector( *( gint->target_A_deltas_pos()),       " target_A_deltas_pos      "); cout << endl;
  cout << endl;

  if ( gint->A_A_deltas_pos()->size() != 0 ) { 

    vector<pair<string,string>> A_A_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->deltas_pos_) )
      A_A_deltas_rngs.push_back(make_pair((*block_rngs_)[ctr.first],(*block_rngs_)[ctr.second])) ;
    
    vector<pair<string,string>> A_A_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->deltas_pos_) ) 
      A_A_deltas_idxs.push_back(make_pair(block_idxs_[ctr.first] , block_idxs_[ctr.second] )) ;


    vector<pair<int,int>> A_A_deltas_pos( gint->deltas_pos_->size() ); 
    vector<pair<int,int>>::iterator aadp_it =  A_A_deltas_pos.begin(); 
    for ( auto gaadp_it = gint->deltas_pos_->begin(); gaadp_it != gint->deltas_pos_->end(); aadp_it++, gaadp_it++ ){
       (*aadp_it).first  =  (*idxs_trans_)[ (*gaadp_it).first  ];
       (*aadp_it).second =  (*idxs_trans_)[ (*gaadp_it).second ];
    }

    vector<pair<string,string>> A_A_deltas_rngs_trans;
    for ( pair<int,int>& ctr :  A_A_deltas_pos )
      A_A_deltas_rngs_trans.push_back(make_pair(std_rngs_[ctr.first] ,std_rngs_[ctr.second])) ;
    
    vector<pair<string,string>> A_A_deltas_idxs_trans;
    for ( pair<int,int>& ctr :  A_A_deltas_pos ) 
      A_A_deltas_idxs_trans.push_back(make_pair((*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] )) ;

    print_pair_vector( A_A_deltas_idxs_trans , "A_A_deltas_idxs_trans " ); 
    print_pair_vector( A_A_deltas_idxs , "   A_A_deltas_idxs " ); 

    assert( A_A_deltas_idxs_trans == A_A_deltas_idxs ); 

  } 

  if ( gint->target_target_deltas_pos()->size() != 0 ) { 

    vector<pair<string,string>> T_T_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) ){
       ctr.first  =  (*idxs_trans_)[ctr.first];
       ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_T_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_T_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_target_deltas_pos()) )
       T_T_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );
  
    print_pair_vector( T_T_deltas_trans_idxs , "T_T_deltas_trans_idxs " ); 
    print_pair_vector( T_T_deltas_idxs , "   T_T_deltas_idxs " ); 

    assert( T_T_deltas_trans_idxs == T_T_deltas_idxs ) ;
  }

  if ( gint->target_A_deltas_pos()->size() != 0 ) { 

 
    vector<pair<string,string>> T_A_deltas_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) ){
      ctr.first  =  (*idxs_trans_)[ctr.first];
      ctr.second =  (*idxs_trans_)[ctr.second];
    }

    vector<pair<string,string>> T_A_deltas_trans_rngs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_trans_rngs.push_back(make_pair( std_rngs_[ctr.first] , std_rngs_[ctr.second] ) );
    
    vector<pair<string,string>> T_A_deltas_trans_idxs;
    for ( pair<int,int>& ctr :  *(gint->target_A_deltas_pos()) )
       T_A_deltas_trans_idxs.push_back(make_pair( (*std_ids_)[ctr.first] , (*std_ids_)[ctr.second] ) );

    print_pair_vector( T_A_deltas_trans_idxs , "T_A_deltas_trans_idxs " ); 
    print_pair_vector( T_A_deltas_idxs , "   T_A_deltas_idxs " ); 

    assert( T_A_deltas_trans_idxs == T_A_deltas_idxs )  ;

 
  } 
  cout << endl << endl;
  return;
}
////////////////////////////////////////////////////////////////////////////
template class GammaGenerator_OrbExcDeriv<double>;
template class GammaGenerator_OrbExcDeriv<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
