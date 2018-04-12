#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Op_General_base::Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                  std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                  std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr< Range_Block_Info> >> all_rxnges ):
                                  idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges),
                                  orig_factor_(factor), num_idxs_(idxs.size()), all_rxnges_(all_rxnges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Op_General_base::Op_General_base constructor"<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  split_rxnges_ = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info> >>();

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Op_General_base::Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                  std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                  std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info> >> split_rxnges ) : 
                                  idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges), orig_factor_(factor), num_idxs_(idxs.size()),
                                  split_rxnges_(split_rxnges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Op_General_base::Op_General_base constructor"<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  all_rxnges_ = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>() ; 

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_Base::TensOp_Base( std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops ) :
                          name_(name), spinfree_(spinfree), Tsymm_("none"), state_dep_(0), 
                          required_blocks_(std::make_shared<std::set<std::string>>()), sub_tensops_(sub_tensops) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
  sort( sub_tensops_.begin(), sub_tensops_.end(), [](shared_ptr<TensOp_Base> t1, shared_ptr<TensOp_Base> t2){ return (bool)( t1->name() < t2->name() );} );
 
  cout << "sub_tensops_ = [ " ; cout.flush();  for ( auto& t : sub_tensops_)    cout << t->name() << " ";  cout << "]" << endl;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is only generated if the indexes of the operator does not range over the same orbitals 
// _not_ have s symmetric ranges, i.e., projection operators. If the operator is (anit-)Hermitian, the results of the map
//  of the map will be obtained from the original all ranges 
//  Only built when needed ( as not necessary for most operators ).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TensOp_Base::TensOp_Base::generate_transformed_ranges( char transformation ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType>::generate_transformed_ranges" <<   endl;

//  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
//  shared_ptr<const vector<bool>> orig_aops = make_shared<const vector<bool>>(aops);

  char trans = tolower(transformation);
  std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info>> trans_ranges_tmp;

  string bob = ""; bob += trans;

  vector<bool> trans_aops = *Op_dense_->aops_ptr_;
  transform_aops_vec( trans, trans_aops );
  print_vector( trans_aops, "trans_aops" ); cout << endl;

  for ( auto& elem : *(Op_dense_->all_rxnges_) ) { 
    vector<string> new_range = elem.first;
    transform_tens_vec( trans, new_range );
    if(satisfies_constraints(new_range)){
      const vector<string> new_range_c = new_range; 
      trans_ranges_tmp.emplace ( new_range_c,  elem.second->get_transformed_block( make_shared<const vector<string>>(new_range_c),
                                                                                   make_shared<const vector<string>>(elem.first),
                                                                                   trans, trans_aops ) ) ; 
      print_vector( new_range_c, bob + " new_range_c" ); cout << endl;
      print_vector( elem.first,  bob + " unique_block" ); cout << endl;
      cout << "got_transformed_block" << endl;
    }
  }
  cout << "got_transformed_block2" << endl;

  if ( trans == 'h' ||  trans == 'i' ) 
    hconj_ranges_ = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info >>> (trans_ranges_tmp);
             
  cout << "got_transformed_block3" << endl;
  return; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Range_Block_Info>
TensOp_Base::TensOp_Base::transform_block_rngs( const vector<char>& rngs, const char op_trans_in ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp_Base::TensOp_Base::transform_block_rngs " << name_ <<  " " << op_trans_in <<  endl;

  char op_trans = tolower(op_trans_in);
 
  pair<double,double> trans_factor;

  vector<char> rngs_tmp = rngs;

  shared_ptr<Range_Block_Info> trans_block; 

  print_vector( rngs_tmp, "rngs_tmp"); cout << endl;

  switch ( op_trans ) { 

    case '0' : 
      trans_block = Op_dense_->all_rxnges_->at( chrvec_to_strvec( rngs_tmp ) );
      break;

    case 'n' : 
      trans_block = Op_dense_->all_rxnges_->at( chrvec_to_strvec( rngs_tmp ) );
      break;

    case 'i' : // inverse
        if ( Op_dense_->hermitian_inverse_ ) {
//          reverse( rngs_tmp.begin(), rngs_tmp.end() );  
          cout << " hello X" << endl;
          vector<string> bob = chrvec_to_strvec( rngs_tmp );
          print_vector( bob , " bob " ); cout << endl;
          cout << " hconj_ranges_->size() = " << hconj_ranges_->size() << endl;
          trans_block = hconj_ranges_->at( chrvec_to_strvec( rngs_tmp ) );
        } else {
          throw logic_error("inverse not defined... aborting!! " ); 
        }
      break;

    case 'h' : // hconj
        reverse( rngs_tmp.begin(), rngs_tmp.end() ); 
        trans_block = hconj_ranges_->at( chrvec_to_strvec( rngs_tmp ) );
      break;

    case 't' : // time reversal
      std::cout << " Why are you time reversing block ranges? Do you mean to time reverse aops_rngs? " << std::endl;
      assert(false);
      break;

    default : 
      string trans_str = "" ; trans_str += op_trans;
      cout << "do not have transformation " << trans_str << " implemented; please check the braket specification in the input file. TensOp_Base::transform_block_rngs" << std::endl;
      assert(false);
      break;
  } 
  print_vector( rngs_tmp , "rngs_tmp after trans " ) ; cout << endl;

  return trans_block;

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                  vector<bool>& aops, DataType orig_factor,
                                  vector<shared_ptr<Transformation>>& symmfuncs,
                                  vector<shared_ptr<Constraint>>& constraints,
                                  string& Tsymm, int state_dep, shared_ptr<map<char, long unsigned int>> range_prime_map ) :
                                  TensOp_Base( name, true ),  symmfuncs_(symmfuncs), constraints_(constraints)   {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp" <<   endl;
          
  CTP_map_ = make_shared< map< string, shared_ptr<CtrTensorPart_Base> >>();
   
  std::vector<int> plus_ops;
  std::vector<int> kill_ops;
  for ( int ii =0 ; ii != aops.size() ; ii++ ) {
    if(aops[ii]) {
      plus_ops.push_back(ii);
    } else {
      kill_ops.push_back(ii);
    }
  }

  auto all_rxnges = generate_rangesX( idxs, idx_ranges, aops ); cout << "left T::generate_rangesX" << endl;

  pair<double,double> orig_factor_tmp =  make_pair(1.0, 1.0);

  Op_dense_ = make_shared<const TensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, orig_factor_tmp, all_rxnges );

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>
TensOp::TensOp<DataType>::generate_rangesX( vector<string>& idxs,  vector<vector<string>>& idx_ranges, vector<bool>& aops ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::TensOp<DataType>::generate_rangesX" <<   endl;

  //set up loop utils
  shared_ptr<vector<int>> no_trans = make_shared<vector<int>>( idxs.size() );
  iota( no_trans->begin(), no_trans->end(), 0 );
  all_ranges_tmp_ =  map< const vector<string>, shared_ptr<Range_Block_Info> >();

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( idxs.size(), 0 );
  shared_ptr<vector<int>> mins = make_shared<vector<int>>( idxs.size(), 0 );
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( idxs.size() );
  pair<double,double> fac_new(1.0,1.0);
  int num_cycles = 1;

  for( int ii = 0; ii != idx_ranges.size(); ii++ ){
    maxs->at(ii) = ( idx_ranges[ii].size()-1 );
    num_cycles *= idx_ranges[ii].size();
  }

  // do all symm test ; loop through ranges , try to transform into another..., if fail,  then add as unique range.
  // generate all possible ranges
  vector<vector<string>> possible_ranges(0);
 do { 
    vector<string> new_block(idx_ranges.size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_block[jj] = idx_ranges[jj][(*fvec)[jj]];

    if(satisfies_constraints(new_block))
      apply_symmetry( new_block, aops ) ;  
     
  } while( fvec_cycle_skipper( fvec, maxs, mins ) );

  return make_shared< const map< const vector<string>, shared_ptr<Range_Block_Info> > > (all_ranges_tmp_);
              
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void  
TensOp::TensOp<DataType>::apply_symmetry( const vector<string>& new_block, std::vector<bool>& aops  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "TensOp::TensOp<DataType>::apply_symmetry" << endl; 

   for ( auto& func : symmfuncs_ ) { 
     auto art_loc = all_ranges_tmp_.find(func->transform( new_block )) ;
     if( art_loc != all_ranges_tmp_.end() ){
       all_ranges_tmp_.emplace( new_block, make_shared<Range_Block_Info>( make_shared<const vector<string>>(new_block) , 
                                                                          art_loc->second->unique_block_,
                                                                          func->idxs_trans(new_block), func->factor(new_block), aops ) ) ;  
       print_vector( new_block , "new_block" ) ; cout << "     "; 
       print_vector( *(art_loc->second->unique_block_), "   unique_block" ) ; cout << "     ";
       print_vector( *(func->idxs_trans(new_block)), "   idxs_trans" ) ; cout << endl;
       break;
     }
   }   

   vector<int> order( aops.size());
   std::iota( order.begin(), order.end(), 0 );
   pair<double,double> new_fac = make_pair( 1.0 , 1.0 );        
   all_ranges_tmp_.emplace( new_block, make_shared<Range_Block_Info>( make_shared<const vector<string>>(new_block) , make_shared<const vector<string>>( new_block ), 
                                                                      make_shared<vector<int>>( order ), new_fac, aops ) );  

   print_vector( new_block,  "new_block"); cout << endl;

   return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
TensOp::TensOp<DataType>::transform_aops( const char op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType>::transform_aops " << endl;

  vector<bool> aops = Op_dense_->aops_; 
  char bob = tolower(op_trans);

  print_vector( aops , "aops original " ) ; cout << endl;

  switch ( bob ) { 

    case 'i' : // inverse
      for_each( aops.begin(), aops.end(), [] ( decltype(aops)::reference aop ) { aop = !aop; } ); 
      break;

    case 'h' : // hconj
      for_each( aops.begin(), aops.end(), [] ( decltype(aops)::reference aop ) { aop = !aop; } ); 
      break;

    case 't' : // time reversal
      break;

    default : 
      string trans_str = "" ; trans_str += op_trans;
      std::cout << "do not have transformation " << trans_str << " implemented; please check the braket specification in the input file. TensOp::transform_aops" << std::endl;
      break;
  } 
  print_vector( aops , "aops after trans " ) ; cout << endl;

  return make_shared<vector<bool>>(aops);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
TensOp_Base::TensOp_Base::transform_aops_rngs( vector<char>& rngs, pair<double,double>& factor, const char op_trans_in ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_Base::TensOp_Base::transform_aops_rngs " << name_ << " " << op_trans_in << endl;
 
  char op_trans = tolower(op_trans_in);

  pair<double,double> trans_factor;

  print_vector( rngs, "aops rngs before trans " ) ; cout << endl;
  switch ( op_trans ) { 

    case 'i' : // inverse
      trans_factor = make_pair( 1.0, -1.0);
      break;

    case 'h' : // hconj
      reverse( rngs.begin(), rngs.end() ); 
      trans_factor = make_pair( 1.0 , -1.0);
      break;

    case 't' : // time reversal
      for_each( rngs.begin(), rngs.end(), [] ( char& rng ) { if ( rng > 'Z' ) { rng -= 32; } else { rng += 32; } } ); 
      trans_factor = make_pair( -1.0, 1.0 ); // determine this from Tsymm of op;
      break;

    default : 
      string trans_str = "" ; trans_str += op_trans;
      std::cout << "do not have transformation " << trans_str << "implemented; please check the braket specification in the input file. TensOp_Base::Transform_aops_rngs" << std::endl;
      break;
  }
  print_vector( rngs, "aops rngs after trans " ) ; cout << endl;

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
MultiTensOp::MultiTensOp<DataType>::transform_aops( const char trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_aops"  << endl;

  throw logic_error( "do not call this yet ; MultiTensOp::MultiTensOp<DataType>::transform_aops , aborting!!"); 

  vector<bool> trans_aops( Op_dense_->num_idxs_, false );
  return make_shared<vector<bool>>(trans_aops);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
TensOp::TensOp<DataType>::transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType>::transform_aops"  << endl;

  throw logic_error( "do not call this yet ;  TensOp::TensOp<DataType>::transform_aops  aborting!!"); 
  vector<bool> trans_aops( Op_dense_->num_idxs_, false );
  return make_shared<vector<bool>>(trans_aops);

}  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<bool>> 
MultiTensOp::MultiTensOp<DataType>::transform_aops( const vector<int>& op_order , const vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_aops " << endl;

  vector<bool> trans_aops( Op_dense_->num_idxs_ );
  vector<bool>::iterator ta_it = trans_aops.begin();
  vector<int>::const_iterator oo_it = op_order.begin(); 
 
  for ( vector<char>::const_iterator ot_it = op_trans.begin(); ot_it != op_trans.end(); ot_it++, oo_it++  ) {

    shared_ptr<vector<bool>> trans_aops_part = sub_tensops_[*oo_it]->transform_aops( *ot_it );
 
    ta_it = move( trans_aops_part->begin(), trans_aops_part->end(), ta_it);

  } 
    
  return make_shared<vector<bool>>(trans_aops);

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<std::vector<char>> 
MultiTensOp::MultiTensOp<DataType>::transform_aops_rngs( shared_ptr<Split_Range_Block_Info> block, pair<double,double>& factor,
                                                         const vector<int>& op_order , const vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_aops_rngs" << endl;

  vector<char> trans_aops_rngs( block->num_idxs_ );
  vector<char>::iterator tar_it = trans_aops_rngs.begin();
  vector<int>::const_iterator oo_it = op_order.begin(); 

  cout << "hello " << endl; 
  for ( vector<char>::const_iterator ot_it = op_trans.begin(); ot_it != op_trans.end(); ot_it++, oo_it++  ) {

   print_vector( *((*block->range_blocks())[*oo_it]->orig_rngs()) , "orig_rngs" ) ; cout << endl;
   vector<char> trans_aops_rngs_part = strvec_to_chrvec( *((*block->range_blocks())[*oo_it]->orig_rngs()) );

   sub_tensops_[*oo_it]->transform_aops_rngs( trans_aops_rngs_part, factor, *ot_it );

   tar_it = move( trans_aops_rngs_part.begin(), trans_aops_rngs_part.end(), tar_it);

  } 

  return make_shared<vector<char>>(trans_aops_rngs);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
std::shared_ptr<Split_Range_Block_Info>
MultiTensOp::MultiTensOp<DataType>::transform_block_rngs( shared_ptr<Split_Range_Block_Info> block,
                                                          shared_ptr<vector<bool>> trans_aops, 
                                                          const vector<int>& op_order , const vector<char>& op_trans ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform_block_rngs" << endl;

  vector<char> trans_block_rngs( block->num_idxs_ );
  vector<char>::iterator tbr_it = trans_block_rngs.begin();
  vector<int>::const_iterator oo_it = op_order.begin(); 

  vector<shared_ptr<Range_Block_Info>> trans_split_block(sub_tensops_.size());
  vector<shared_ptr<Range_Block_Info>>::iterator tsb_it = trans_split_block.begin();

  vector<int> rngs_trans_merged( block->num_idxs_);

  cout << "hello " << endl; 
  for ( vector<char>::const_iterator ot_it = op_trans.begin(); ot_it != op_trans.end(); ot_it++, oo_it++ , tsb_it++ ) {

    print_vector( *((*block->range_blocks())[*oo_it]->orig_rngs()) , "orig_rngs" ) ; cout << endl;
    vector<char> trans_block_rngs_part = strvec_to_chrvec( *((*block->range_blocks())[*oo_it]->orig_rngs()) );
    
    *tsb_it = sub_tensops_[*oo_it]->transform_block_rngs( trans_block_rngs_part, *ot_it );
    
    tbr_it = move( trans_block_rngs_part.begin(), trans_block_rngs_part.end(), tbr_it);
  } 
  SRBI_Helper helper( trans_split_block );
  auto srbi = make_shared<Split_Range_Block_Info>( make_shared<const vector<bool>>(*trans_aops) , helper );
 
  print_vector( trans_block_rngs, "trans_block_rngs" ) ; cout << endl; 
  
  return srbi; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
bool TensOp::TensOp<DataType>::satisfies_constraints( vector<string>& ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (auto& cstr : constraints_) 
    if(!(cstr->apply_constraint( ranges )) )
      return false;
  return true;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and ranges; note this is done by TensOp not TensOp gen as the contractions needed will be 
//dependent on the sparsity associated with the relevant state
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::generate_uncontracted_ctps() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType> generate_uncontracted_ctps" /* << name_*/ << endl;

  //TODO fix this hack; only call get_ctrs_tens_ranges after a shared pointer has been created elsewhere, and get_ctrs_tens_ranges
  //     has to be called before anything can happen to the TensOp, so this is basically safe, just totally illogical.
  sub_tensops_ = vector<shared_ptr<TensOp_Base>>(1 , this->shared_from_this());

  //puts uncontracted ranges into map 
  shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( Op_dense_->idxs_); // change this to not a const here or to a const in ctp;
  shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector< pair<int,int>>>(0);
  for (auto rng_it = all_rxnges()->begin(); rng_it != all_rxnges()->end(); rng_it++) {
    shared_ptr<vector<pair<int,int>>>  ReIm_factors = make_shared< vector<pair<int,int>>>(1, rng_it->second->factors()); 
    shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
 
    shared_ptr<CtrTensorPart<DataType>>  CTP = make_shared< CtrTensorPart<DataType> >( full_idxs, full_ranges, noctrs, ReIm_factors ); 
    cout << "CTP->name() = " << CTP->name() << " is going into CTP_map" <<  endl;
    CTP_map_->emplace(CTP->name(), CTP); //maybe should be addded in with ctr_idxs
  }
  return;
} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
MultiTensOp::MultiTensOp<DataType>::MultiTensOp( std::string name, bool spinfree,
                                                 std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops,
                                                 shared_ptr<map< char , long unsigned int>> range_prime_map ):
                                                 TensOp_Base( name, spinfree, sub_tensops ),
                                                 num_tensors_(sub_tensops_.size()) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::MultiTensOp" << endl;
    
  vector<string> idxs;
  vector<vector<string>> idx_ranges;
  vector<bool> aops;
 
  vector<int> cmlsizevec(num_tensors_);
  vector<int> plus_ops;
  vector<int> kill_ops;
 
  int num_idxs = 0;
  for ( int ii = 0;  ii != sub_tensops_.size() ; ii++ ) {
  
    cmlsizevec[ii]  = num_idxs;
    idx_ranges.insert(  idx_ranges.end(), sub_tensops_[ii]->idx_ranges()->begin(), sub_tensops_[ii]->idx_ranges()->end() );
    idxs.insert( idxs.end(),sub_tensops_[ii]->idxs()->begin(), sub_tensops_[ii]->idxs()->end() );
    aops.insert( aops.end(),sub_tensops_[ii]->aops()->begin(), sub_tensops_[ii]->aops()->end() );

    for ( int elem : *(sub_tensops_[ii]->plus_ops()) )
      plus_ops.push_back(elem + cmlsizevec[ii]);
  
    for ( int elem : *(sub_tensops_[ii]->kill_ops()) )
      kill_ops.push_back(elem + cmlsizevec[ii]);
  
    num_idxs += sub_tensops_[ii]->num_idxs();
  } 
  shared_ptr<const map < const vector<string> , shared_ptr<Split_Range_Block_Info > > >  all_rxnges_ptr = generate_rangesX( idxs, aops, cmlsizevec );
 
  Op_dense_ = make_shared<const MultiTensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, make_pair(1.0,1.0), cmlsizevec, all_rxnges_ptr ); 
  CTP_map_  = make_shared< map< string, shared_ptr<CtrTensorPart_Base> >>();

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr< const map< const vector<string>, shared_ptr<Split_Range_Block_Info> >>
MultiTensOp::MultiTensOp<DataType>::generate_rangesX( vector<string>& idxs, vector<bool>& aops, vector<int>& cmlsizevec  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_rangesX()" << endl;
 
  int num_idxs = idxs.size();
  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
  vector< map< const  vector<string>, shared_ptr< Range_Block_Info >>::const_iterator> rng_maps(num_tensors_);  
  vector<int> posvec( num_tensors_, 0 ); 

  map< const vector<string>, shared_ptr<Split_Range_Block_Info> > split_rxnges;

  if ( num_tensors_ > 1 ) { 
    shared_ptr<vector<int>> forvec = make_shared<vector<int>>(num_tensors_, 0 ); 
    shared_ptr<vector<int>> old_forvec = make_shared<vector<int>>(*forvec);
 
    shared_ptr<vector<int>> mins   = make_shared<vector<int>>(num_tensors_, 0 );  
    shared_ptr<vector<int>> maxs   = make_shared<vector<int>>(num_tensors_ );  

    for ( int ii = 0 ; ii != sub_tensops_.size() ; ii++ ){ 
      maxs->at(ii) = sub_tensops_[ii]->all_rxnges()->size()-1;
      rng_maps[ii] = sub_tensops_[ii]->all_rxnges()->begin() ; 
    }
 
    do {
 
      for ( int ii = 0 ; ii != forvec->size() ; ii++ ) {
        if ( (*old_forvec)[ii] != (*forvec)[ii] ) {
          if ( (*forvec)[ii] == 0 ) {
            rng_maps[ii] = sub_tensops_[ii]->all_rxnges()->begin() ;
          } else {
            rng_maps[ii]++;
          }
        }
       }
     
      old_forvec = make_shared<vector<int>>(*forvec);
      shared_ptr< vector <shared_ptr<Range_Block_Info >>> split_block = make_shared< vector <shared_ptr<Range_Block_Info >>>( num_tensors_ );
      vector<shared_ptr<Range_Block_Info>>::iterator sb_it = split_block->begin();
      for (auto  rm_it = rng_maps.begin(); rm_it != rng_maps.end(); rm_it++, sb_it++ )
        *sb_it = (*rm_it)->second;

      //TODO Must obtain from constraint functions 
      shared_ptr<Split_Range_Block_Info> srbi;
      {
        SRBI_Helper helper(*split_block);
        srbi = make_shared<Split_Range_Block_Info>( make_shared<const vector<bool>>(aops), helper );
      }
      split_rxnges.emplace( *(srbi->orig_rngs()), srbi ) ;

    } while( fvec_cycle_skipper( forvec, maxs, mins ) );

    
  } else { 

    //TODO add constructor so can just have single tensor, need not be vector
    for ( auto elem : *(sub_tensops_[0]->all_rxnges()) )  {
      SRBI_Helper helper( make_shared<vector<shared_ptr<Range_Block_Info >>>(1, elem.second));
      shared_ptr<Split_Range_Block_Info> srbi = make_shared<Split_Range_Block_Info>( Op_dense_->aops(), helper );
      split_rxnges.emplace( elem.first, srbi ) ;
    }
  }
   
  cout << "split_rxnges.size() = " << split_rxnges.size() << endl;
  return make_shared< const map < const vector<string> , shared_ptr<Split_Range_Block_Info > > >( split_rxnges );
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::generate_uncontracted_ctps() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_uncontracted_ctps " <<  endl;

  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);
  for (auto rng_it = Op_dense_->split_rxnges()->begin(); rng_it != Op_dense_->split_rxnges()->end(); rng_it++) 
    enter_cmtps_into_map(*noctrs, rng_it->second->factors(), rng_it->first );
  
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::enter_cmtps_into_map(pint_vec ctr_pos_list, pair<int,int> ReIm_factors, const vector<string>& id_ranges ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // cout << "MultiTensOp::enter_into_CTP_map" << endl;

  shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>> CTP_vec = make_shared< vector< shared_ptr<CtrTensorPart_Base> >> (num_tensors_); 
  vector<pair<pair<int,int>,pair<int,int>>> diffT_ctrs_pos(0);
  vector<vector<pair<int,int>>> sameT_ctrs_pos( num_tensors_,  pint_vec(0));
  shared_ptr<vector<pair<int,int>>> no_ctrs =  make_shared<vector<pair<int,int>>>(0);
  shared_ptr<vector<pair<int,int>>> ReIm_factor_vec = make_shared<vector<pair<int,int>>>(1, ReIm_factors ) ; 

  //seperate contractions into those on the same tensor, and those between different tensors 
  // TODO tidy this up, e.g., use lambda to get cross pos
  for ( pair<int,int> ctr_pos : ctr_pos_list ) {
 
    pair<int,int> ctr1;
    pair<int,int> ctr2;

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.first -= sub_tensops_[ii]->num_idxs()) < 0 ){
         ctr1 = make_pair(ii, ctr_pos.first + sub_tensops_[ii]->num_idxs());
        break;
      }

    for ( int ii = 0; ii != num_tensors_ ; ii++ )  
      if ( (ctr_pos.second -= sub_tensops_[ii]->num_idxs()) < 0 ){
         ctr2 = make_pair(ii, ctr_pos.second + sub_tensops_[ii]->num_idxs());
        break;
      } 

    if (ctr1.first == ctr2.first) {
      sameT_ctrs_pos.at(ctr1.first).push_back(make_pair(ctr1.second, ctr2.second));  
    } else {
      diffT_ctrs_pos.push_back(make_pair(ctr1,ctr2));
    }
  } 

  //get_ranges for individual tensors
  for (int ii = 0 ; ii !=num_tensors_; ii++ ){ 
    shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(id_ranges.begin() + Op_dense_->cmlsizevec(ii), id_ranges.begin()+cmlsizevec(ii)+sub_tensops_[ii]->num_idxs());
    shared_ptr<vector<string>>  TS_idxs = make_shared<vector<string>>( *(sub_tensops_[ii]->idxs()) );
   
    if( sameT_ctrs_pos[ii].size() != 0 ) {
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, make_shared<vector<pair<int,int>>>(sameT_ctrs_pos[ii]), ReIm_factor_vec ) ; 
    } else { 
      CTP_vec->at(ii) = make_shared< CtrTensorPart<DataType> >( TS_idxs, TS_id_ranges, no_ctrs, ReIm_factor_vec ); 
    }
    CTP_map_->emplace(CTP_vec->at(ii)->name(), CTP_vec->at(ii)); 

  }
  if ( CTP_vec->size() < 3  ){ 
    auto new_cmtp = make_shared<CtrMultiTensorPart<DataType>>(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos));
    CTP_map_->emplace(new_cmtp->name(), new_cmtp );
  } else {
    get_cmtp(CTP_vec, make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(diffT_ctrs_pos) ); 
  } 

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If the ctp_vec has more than two elements then the multitensor is built up recursively from other multitensors, this is
// to simplify the process of contracting the tensors.  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::get_cmtp( shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>>  ctp_vec,  // ctp : contracted_tensor_parts
                                                   shared_ptr<vector<pair<pair<int,int>, pair<int,int>>>> ccp_vec  ){ // ccp : cross_ctrs_pos
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // cout << "MultiTensOp::MultiTensOp<DataType>::get_cmtp " << endl;
  
  auto ctp_vec_buff =  make_shared<vector<shared_ptr<CtrTensorPart_Base>>>(*ctp_vec);  
  auto ccp_vec_buff =  make_shared<vector<pair<pair<int,int>, pair<int,int>>>>(*ccp_vec);
  auto cmtp_orig_ctp_order = make_shared<CtrMultiTensorPart<DataType>>(ctp_vec_buff, ccp_vec_buff);

  vector<pair<pair<int,int>, pair<int,int>>>::iterator ccp_it = ccp_vec->begin();
  int counter = 0 ;
  if ( ccp_vec->size() != 0 ) { 
    do {
      int inside_counter = 0 ;
      int ta_pos = ccp_it->first.first < ccp_it->second.first ? ccp_it->first.first : ccp_it->second.first ;
      int tb_pos = ccp_it->first.first > ccp_it->second.first ? ccp_it->first.first : ccp_it->second.first ;

      auto ctp_vec_tatb = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( vector<shared_ptr<CtrTensorPart_Base>> { ctp_vec->at(ta_pos), ctp_vec->at(tb_pos) } );
      auto ccp_vec_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);
      auto ccp_vec_merged_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);

      for ( auto ccp_ab_it = ccp_vec->begin(); ccp_ab_it != ccp_vec->end();  ccp_ab_it++ ) {
        if (ccp_ab_it->first.first == ta_pos && ccp_ab_it->second.first == tb_pos) {
          ccp_vec_tatb->push_back( make_pair( make_pair( 0, ccp_ab_it->first.second), make_pair(1, ccp_ab_it->second.second) ));
        } else if (ccp_ab_it->first.first == tb_pos && ccp_ab_it->second.first == ta_pos){
          ccp_vec_tatb->push_back( make_pair( make_pair( 0, ccp_ab_it->second.second), make_pair(1, ccp_ab_it->first.second) ));
        } else {
          ccp_vec_merged_tatb->push_back(*ccp_ab_it);
        }
      }

      shared_ptr<CtrMultiTensorPart<DataType>> cmtp_tatb = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec_tatb, ccp_vec_tatb );
      CTP_map_->emplace( cmtp_tatb->name(), cmtp_tatb );

      shift_ccp_and_ctp_vecs( cmtp_tatb, ta_pos, tb_pos, ctp_vec, ccp_vec_merged_tatb );

      ccp_vec = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(*ccp_vec_merged_tatb);
      ccp_it = ccp_vec->begin();

    } while( ccp_vec->size() != 0 ) ;

  } else { 
    while ( ctp_vec->size() > 2 ) { 
      auto ctp_vec_tatb = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( vector<shared_ptr<CtrTensorPart_Base>> { *(ctp_vec->end()-2), *(ctp_vec->end()-1) } );
      auto ccp_vec_tatb = make_shared<vector<pair<pair<int,int>, pair<int,int>> >>(0);
      auto cmtp_tatb = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec_tatb, ccp_vec_tatb );
      CTP_map_->emplace( cmtp_tatb->name(), cmtp_tatb );
      ctp_vec->pop_back(); 
      ctp_vec->back() = cmtp_tatb;
    }
  }   

  shared_ptr<CtrMultiTensorPart<DataType>> cmtp_new_ctp_order = make_shared<CtrMultiTensorPart<DataType>>( ctp_vec, ccp_vec );
  std::shared_ptr<std::vector<int>> new_order_to_orig_order =  get_pattern_match_order( cmtp_orig_ctp_order->unc_idxs(), cmtp_new_ctp_order->unc_idxs() ) ; 

  cmtp_orig_ctp_order->use_new_order_compute_list( new_order_to_orig_order, cmtp_new_ctp_order->name() );

  CTP_map_->emplace( cmtp_orig_ctp_order->name(), cmtp_orig_ctp_order ) ;
  CTP_map_->emplace( cmtp_new_ctp_order->name(), cmtp_new_ctp_order );

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Shifts the indexes in the cross_ctrs_vec (necessary as ctps at ta and tb have been merged and put to the end or ctp_vec)
// Note that both ctp vec and ccp vec will be altered in this routine
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
MultiTensOp::MultiTensOp<DataType>::shift_ccp_and_ctp_vecs( shared_ptr<CtrMultiTensorPart<DataType>>& tatb_cmtp,
                                                            int ta, int tb, shared_ptr<vector<shared_ptr<CtrTensorPart_Base>>>& ctp_vec,
                                                            shared_ptr<vector<pair<pair<int,int>, pair<int,int>> >>& ccp_vec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "MultiTensOp::MultiTensOp<DataType>::shift_ccp_ctp_vecs" << endl;

  auto new_ctp_vec = make_shared<vector<shared_ptr<CtrTensorPart_Base>>>( ctp_vec->size()-1 );
  auto new_ctp_vec_it = new_ctp_vec->begin();
  map<int,int> shifted_ctp_pos_map;
  int tb_id_shift = ctp_vec->at(ta)->size();

  int new_ctp_pos = 0;
  for (int ii = 0 ; ii != ctp_vec->size() ; ii++ ){
    if ( ii != ta && ii != tb ){
      *new_ctp_vec_it++ = ctp_vec->at(ii);
       shifted_ctp_pos_map.emplace(ii,new_ctp_pos++);
    }
  }
  new_ctp_vec->back() = tatb_cmtp;
  ctp_vec = new_ctp_vec;

  shifted_ctp_pos_map.emplace(ta, new_ctp_vec->size()-1);
  shifted_ctp_pos_map.emplace(tb, new_ctp_vec->size()-1);

  for ( pair<pair<int,int> ,pair<int,int>>& ccp : *ccp_vec ) {
    if (ccp.first.first == tb ) {
      ccp.first.second += tb_id_shift;
    } else if (ccp.second.first == tb ) {
      ccp.second.second += tb_id_shift;
    }
    ccp.first.first = shifted_ctp_pos_map.at(ccp.first.first);
    ccp.second.first = shifted_ctp_pos_map.at(ccp.second.first);
  }

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp::TensOp<double>;
template class TensOp::TensOp<complex<double>>;
template class MultiTensOp::MultiTensOp<double>;
template class MultiTensOp::MultiTensOp<complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
