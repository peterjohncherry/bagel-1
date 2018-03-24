#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/tensop.h>
#include <src/prop/proptool/algebraic_manipulator/algebra_utils.h>

using namespace std;
using namespace WickUtils;
using namespace Algebra_Utils;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Op_General_base::Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                  std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                  std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                  std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr< Range_Block_Info> >> all_ranges ):
                                  idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges), orig_factor_(factor), num_idxs_(idxs.size()),
                                  unique_range_blocks_(*unique_range_blocks), all_ranges_(all_ranges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Op_General_base::Op_General_base constructor"<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  unique_range_blocks_ptr_ = make_shared<const vector< shared_ptr< const vector<string>>>>(unique_range_blocks_);
  
  split_ranges_ = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info> >>();

  all_ranges_trans_ = make_shared<map< char, shared_ptr<const map< const vector<string>,  shared_ptr<Range_Block_Info>>> >>();

  split_ranges_trans_ = make_shared< map< pair<vector<char>, vector<int>>, shared_ptr<const map< const vector<string>,  shared_ptr<Split_Range_Block_Info>>> >>();

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Op_General_base::Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                  std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                  std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks ,
                                  std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info> >> split_ranges ) : 
                                  idxs_(idxs),  aops_(aops), plus_ops_(plus_ops), kill_ops_(kill_ops), idx_ranges_(idx_ranges), orig_factor_(factor), num_idxs_(idxs.size()),
                                  unique_range_blocks_(*unique_range_blocks), split_ranges_(split_ranges) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Op_General_base::Op_General_base constructor"<< endl;
           
  idxs_ptr_ =  make_shared<const vector<string>>(idxs_);
  aops_ptr_ =  make_shared<const vector<bool>>(aops_);
  idx_ranges_ptr_ =  make_shared<const vector<vector<string>>>(idx_ranges_);

  plus_ops_ptr_ = make_shared<const vector<int>>(plus_ops);
  kill_ops_ptr_ = make_shared<const vector<int>>(kill_ops);

  unique_range_blocks_ptr_ = make_shared<const vector< shared_ptr< const vector<string>>>>(unique_range_blocks_);

  all_ranges_ = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>() ; 

  all_ranges_trans_ = make_shared<map< char, shared_ptr<const map< const vector<string>,  shared_ptr<Range_Block_Info>>> >>();

  split_ranges_trans_ = make_shared< map< pair<vector<char>, vector<int>>, shared_ptr<const map< const vector<string>,  shared_ptr<Split_Range_Block_Info>>> >>();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_General::TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info >>> split_ranges) : 
                                Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, unique_range_blocks, split_ranges) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_General::TensOp_General "<< endl;

   // no split ranges
   // Does there even need to be anything here? 
   // Yes, but can leave it for now; must template general_ops and Range_Block_Info; hermitian conjugation will mess things around,
   // but want all the basic functions in base class. Presumably the one for returning factor, unique block and all_ranges will
   // be virtual there, and only defined properly here... Will need range_block base as well I guess.

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TensOp_General::TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                                std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info >>> all_ranges) : 
                                Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, unique_range_blocks, all_ranges ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_General::TensOp_General "<< endl;
  // has split ranges
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MultiTensOp_General::MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                                          std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                                          std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                                          std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr< Split_Range_Block_Info>>> split_ranges ) :
                                          Op_General_base::Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, unique_range_blocks, split_ranges ), 
                                          cmlsizevec_(cmlsizevec) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "MultiTensOp_General::MultiTensOp_General "<< endl;
  //  all_ranges_ =  make_shared<const map< const vector<string>, shared_ptr<Split_Range_Block_Info >>>(all_ranges_);
  cmlsizevec_ptr_ =  make_shared<const vector<int>>(cmlsizevec_);
  split_ranges_ = split_ranges; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
TensOp::TensOp<DataType>::TensOp( string name, vector<string>& idxs, vector<vector<string>>& idx_ranges,
                                  vector<bool>& aops, DataType orig_factor,
                                  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>), int, int > >& symmfuncs, 
                                  vector<bool(*)(shared_ptr<vector<string>>) >& constraints,
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

  auto test_var = generate_ranges( idxs, idx_ranges, aops, range_prime_map );

  pair<double,double> orig_factor_tmp =  make_pair(1.0, 1.0);

  Op_dense_ = make_shared<const TensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, orig_factor_tmp, get<1>(test_var), get<0>(test_var) );

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO this is too slow; just build ranges using symmetry, don't check.
template<typename DataType>
tuple< std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>,
       std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>> >
TensOp::TensOp<DataType>::generate_ranges( vector<string>& idxs,  vector<vector<string>>& idx_ranges, vector<bool>& aops,
                                           shared_ptr<map< char , long unsigned int >> range_prime_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::generate_ranges" <<   endl;

  //set up loop utils

  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
  int num_idxs = orig_idxs->size();

  vector< shared_ptr< const vector<string>> > unique_range_blocks;
  map< const vector<string>, shared_ptr<Range_Block_Info> > all_ranges;
  
  auto apply_symmetry = [ &all_ranges, &orig_idxs, &aops, &range_prime_map ]( const vector<string>& ranges_1,  const vector<string>& ranges_2 ) {
     const bool a_unique_range_block = true; //TODO put symmetry back in
     if (a_unique_range_block) {
       const pair<int,int> fac(1,1);
       const pair<double,double> fac_new(1.0,1.0);
       shared_ptr<const vector<string>> unique_ranges = make_shared<const vector<string>>(ranges_2); 
       shared_ptr<const vector<string>> trans_idxs = orig_idxs;
       bool survives = WickUtils::RangeCheck( ranges_2, aops );
       shared_ptr<const vector<bool>> test = make_shared<const vector<bool>>(aops); 
       all_ranges.emplace(ranges_2, make_shared< Range_Block_Info >( a_unique_range_block, survives, fac_new, unique_ranges, unique_ranges, orig_idxs, trans_idxs, test, range_prime_map ) );
     }
     return false;
  };

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( num_idxs, 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( num_idxs );
  bool is_unique = true;
  pair<int,int> fac(1,1);
  pair<double,double> fac_new(1.0,1.0);
  int  num_cycles = 1;

  for( int ii = 0; ii != idx_ranges.size(); ii++ ){
    maxs->at(ii) = ( idx_ranges[ii].size()-1 );
    num_cycles *= idx_ranges[ii].size();
  }

  //generate all possible ranges
  vector<vector<string>> possible_ranges(0);
  for (int ii = 0 ; ii != num_cycles ; ii++){
    vector<string> new_range(idx_ranges.size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range[jj] = idx_ranges[jj][fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges.push_back(new_range);
    fvec_cycle( fvec, maxs );
  } 

  // Initialize range maps and add in first range to start things off.
  const vector<string> init_range = possible_ranges[0];
  shared_ptr<const vector<string>> orig_range = make_shared<const vector<string>>(possible_ranges[0]);
  shared_ptr<const vector<string>> trans_range = orig_range;
  shared_ptr<const vector<string>> trans_idxs = orig_idxs;
   
  bool survives = WickUtils::RangeCheck( *orig_range, aops );
  shared_ptr<const vector<bool>> test = make_shared<const vector<bool>>(aops); 
  all_ranges.emplace(init_range, make_shared< Range_Block_Info >( is_unique, survives, fac_new, orig_range, trans_range, orig_idxs, trans_idxs, test, range_prime_map ) );
 
  unique_range_blocks = vector< shared_ptr< const vector<string>>>(1, trans_range );
  //Apply symmetry operations to remove unnecessary ranges 
  for (int ii = 1 ; ii!=possible_ranges.size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
      if( apply_symmetry(possible_ranges[kk], possible_ranges[ii]) ){ //note that this sets the elements of all_ranges_
        break;
      }
      if(kk == ii-1){
        unique_range_blocks.push_back(make_shared<const vector<string>>(possible_ranges[ii]));
        break;
      }
    } 
  }

  auto all_ranges_ptr = make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info >>> (all_ranges);
  auto unique_ranges_ptr = make_shared<std::vector< std::shared_ptr< const std::vector<std::string>> >>(unique_range_blocks);                                                                                                                
  return tie ( all_ranges_ptr, unique_ranges_ptr );
              
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::transform( char transformation,  shared_ptr<map< char , long unsigned int >> range_prime_map ) {  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is not done in range_block_info as the symmetry functions are stored in here; note that the transformations performed by the symmetry functions
// pertain to the transformations of the elements of the molecular orbital integral tensor. They *do not* transform the creation and annihilation 
// operators. The operator transformations performed in this function *do* transform the creation and annihilation operators. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Range_Block_Info::transform( bool is_unique, bool survives, std::pair<double,double> factors " << endl; 
  //set up loop utils

  vector<string> trans_orig_idxs = *idxs(); 
  vector<bool> trans_orig_aops = *aops(); 
  transform_tens_vec( transformation, trans_orig_idxs ); 
  transform_tens_vec( transformation, trans_orig_aops ); 
 
  int num_idxs = trans_orig_idxs.size();
  vector< shared_ptr< const vector<string>> > unique_range_blocks;
  map< const vector<string>, shared_ptr<Range_Block_Info> > all_ranges_tmp;
  
  auto apply_symmetry = [&all_ranges_tmp,  &trans_orig_idxs, &trans_orig_aops, &range_prime_map, this ]( const vector<string>& ranges_1,  const vector<string>& ranges_2 ) {
     bool a_unique_range_block = true; //TODO put symmetry back in, should compare to original all ranges....
     if (a_unique_range_block) {
       pair<double,double> block_fac(1.0,1.0); // TODO set from transform
       shared_ptr<const vector<string>> unique_ranges_trans = make_shared<const vector<string>>(ranges_2); 
       shared_ptr<const vector<string>> trans_orig_idxs_ptr = make_shared<const vector<string>>(trans_orig_idxs);
       bool survives = WickUtils::RangeCheck( ranges_2, trans_orig_aops );
       shared_ptr<const vector<bool>> trans_orig_aops_ptr = make_shared<const vector<bool>>( trans_orig_aops ); 
       all_ranges_tmp.emplace(ranges_2, make_shared< Range_Block_Info >( a_unique_range_block, survives, block_fac, unique_ranges_trans, unique_ranges_trans, trans_orig_idxs_ptr, trans_orig_idxs_ptr, trans_orig_aops_ptr, range_prime_map ) );
     }
     return false;
  };

  shared_ptr<vector<int>> fvec = make_shared<vector<int>>( num_idxs, 0); 
  shared_ptr<vector<int>> mins = make_shared<vector<int>>( num_idxs, 0); 
  shared_ptr<vector<int>> maxs = make_shared<vector<int>>( num_idxs );
  int  num_cycles = 1;
  vector<vector<string>> trans_idx_ranges = *idx_ranges();
  for( int ii = 0; ii != trans_idx_ranges.size(); ii++ ){
    maxs->at(ii) = ( trans_idx_ranges[ii].size()-1 );
    num_cycles *= trans_idx_ranges[ii].size();
  }

  //generate all possible ranges
  vector<vector<string>> possible_ranges(0);
  for (int ii = 0 ; ii != num_cycles ; ii++){
    vector<string> new_range(trans_idx_ranges.size());
    for (auto jj = 0 ; jj != fvec->size() ; jj++ )
       new_range[jj] = trans_idx_ranges[jj][fvec->at(jj)];
    if(satisfies_constraints(new_range))
      possible_ranges.push_back(new_range);
    fvec_cycle_skipper( fvec, maxs, mins );
  } 

  //Apply symmetry operations to remove unnecessary ranges 
  for (int ii = 0 ; ii!=possible_ranges.size(); ii++) {
    for (int kk = 0 ; kk != ii; kk++) {
      if( apply_symmetry(possible_ranges[kk], possible_ranges[ii]) ){ //note that this sets the elements of all_ranges_trans_
        break;
      }
      if(kk == ii-1){
        unique_range_blocks.push_back(make_shared<const vector<string>>(possible_ranges[ii]));
        break;
      }
    } 
  }


  Op_dense_->all_ranges_trans_->emplace( transformation, make_shared<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info >>> (all_ranges_tmp) );
  return;
  
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
bool TensOp::TensOp<DataType>::satisfies_constraints( vector<string>& ranges ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  shared_ptr<vector<string>> ranges_ptr = make_shared<vector<string>>(ranges);
  for (auto OK : constraints_) 
    if( !OK(ranges_ptr))
      return false;
  return true;
}; 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and ranges; note this is done by TensOp not TensOp gen as the contractions needed will be 
//dependent on the sparsity associated with the relevant state
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp::TensOp<DataType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef DBG_TensOp 
  cout << "TensOp get_ctrs_tens_ranges" << endl;
#endif 
//////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp::TensOp<DataType> get_ctrs_tens_ranges" /* << name_*/ << endl;

  //TODO fix this hack; only call get_ctrs_tens_ranges after a shared pointer has been created elsewhere, and get_ctrs_tens_ranges
  //     has to be called before anything can happen to the TensOp, so this is basically safe, just totally illogical.
  sub_tensops_ = vector<shared_ptr<TensOp_Base>>(1 , this->shared_from_this());

  //puts uncontracted ranges into map 
  shared_ptr<vector<pair<int,int>>>  noctrs = make_shared<vector< pair<int,int>>>(0);
  for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
    shared_ptr<vector<pair<int,int>>>  ReIm_factors = make_shared< vector<pair<int,int>>>(1, rng_it->second->factors()); 
    shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
    shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );
 
    shared_ptr<CtrTensorPart<DataType>>  CTP = make_shared< CtrTensorPart<DataType> >( full_idxs, full_ranges, noctrs, ReIm_factors ); 
    cout << "CTP->name() = " << CTP->name() << " is going into CTP_map" <<  endl;
    CTP_map_->emplace(CTP->name(), CTP); //maybe should be addded in with ctr_idxs
  }
 
  //puts_contracted ranges into map
  for ( int nctrs = 1 ; nctrs != (num_idxs()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>>  ctr_lists = get_unique_pairs(plus_ops(), kill_ops(), nctrs);
    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = all_ranges()->begin(); rng_it != all_ranges()->end(); rng_it++) {
        bool valid =true;
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[ctr_vec->at(ii).first] != rng_it->first[ctr_vec->at(ii).second]){
            valid = false;
            break;
          }
        }

        if (valid) {
          shared_ptr<vector<pair<int,int>>> ReIm_factors = make_shared<vector<pair<int,int>>>(1, rng_it->second->factors()); 
          shared_ptr<vector<string>> full_ranges = make_shared<vector<string>>(rng_it->first);
          shared_ptr<vector<string>> full_idxs   = make_shared<vector<string>>( *this->idxs() );

          shared_ptr<CtrTensorPart<DataType>> CTP = make_shared<CtrTensorPart<DataType>>( full_idxs, full_ranges, ctr_vec, ReIm_factors ); 
          cout << "CTP->name() = " << CTP->name() << "is going into CTP_map " << endl;
          CTP_map_->emplace(CTP->name(), CTP); 

        }
      }
    }
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Only build the CTP for the unique ranges; these are the only ones needed. If gamma generator
// asks for any other ones, then you have a bug.
///////////////////////////.///////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DType>
void TensOp::TensOp<DType>::get_ctp_idxs_ranges( shared_ptr<vector<pair<int,int>>> ctrs_pos, shared_ptr<Range_Block_Info> block_info ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<bool> get_unc(block_info->orig_idxs()->size(), true);
  for (int ii =0; ii != ctrs_pos->size() ; ii++){
    get_unc[ctrs_pos->at(ii).first] = false;
    get_unc[ctrs_pos->at(ii).second] = false;
  }

  int num_unc_ids =  get_unc.size() - ctrs_pos->size()*2;

  vector<string> unc_id_ranges( num_unc_ids );
  vector<string> unc_idxs( num_unc_ids );
  vector<int> unc_pos( num_unc_ids );
  map<int,int> unc_rel_pos;

  int jj = 0;
  for ( int ii = 0 ; ii !=get_unc.size() ; ii++ ) {
    if (get_unc[ii]){
      unc_id_ranges[jj] = block_info->orig_block()->at(ii);
      unc_idxs[jj]      = block_info->orig_idxs()->at(ii);
      unc_pos[jj]       = ii;
      unc_rel_pos.emplace(ii, jj);
      jj++;
    }
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
  shared_ptr<const map < const vector<string> , shared_ptr<Split_Range_Block_Info > > >  all_ranges_ptr = generate_ranges( idxs, aops, cmlsizevec, range_prime_map );

  auto  unique_range_blocks = make_shared<std::vector< std::shared_ptr<const std::vector<std::string>>>>(0);
  for ( auto all_range_map_it : *all_ranges_ptr ) 
    unique_range_blocks->push_back( make_shared<const vector<string>>(all_range_map_it.first) ); 
  
  Op_dense_ = make_shared<const MultiTensOp_General>( idxs, aops, plus_ops, kill_ops, idx_ranges, make_pair(1.0,1.0), cmlsizevec, unique_range_blocks, all_ranges_ptr ); 
  CTP_map_  = make_shared< map< string, shared_ptr<CtrTensorPart_Base> >>();

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
shared_ptr< const map< const vector<string>, shared_ptr<Split_Range_Block_Info> >>
MultiTensOp::MultiTensOp<DataType>::generate_ranges( vector<string>& idxs, vector<bool>& aops, vector<int>& cmlsizevec, 
                                                     shared_ptr<map<char, long unsigned int>> range_prime_map ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::generate_ranges()" << endl;

  int num_idxs = idxs.size();
  shared_ptr<const vector<string>> orig_idxs = make_shared<const vector<string>>(idxs);
  vector< map< const  vector<string>, shared_ptr< Range_Block_Info >>::const_iterator> rng_maps(num_tensors_);  
  vector<int> posvec( num_tensors_, 0 ); 

  map< const vector<string>, shared_ptr<Split_Range_Block_Info> > all_ranges;

  // Such an ugly loop that this should be fixed when switch ranges to numbers
  // The reason for it is that maps can't map->begin()+ii , so must iterate
  // This is why you had the ugly thing before; it was much faster than this stupid thing with all it's checking.
  // Maybe add some comments in future so you can understand things...
  if ( num_tensors_ > 1 ) { 
    shared_ptr<vector<int>> forvec = make_shared<vector<int>>(num_tensors_, 0 ); 
    shared_ptr<vector<int>> old_forvec = make_shared<vector<int>>(*forvec);

    shared_ptr<vector<int>> mins   = make_shared<vector<int>>(num_tensors_, 0 );  
    shared_ptr<vector<int>> maxs   = make_shared<vector<int>>(num_tensors_, 0 );  
    for ( int ii = 0 ; ii != sub_tensops_.size() ; ii++ ){ 
      maxs->at(ii) = sub_tensops_[ii]->all_ranges()->size()-1;
      rng_maps[ii] = sub_tensops_[ii]->all_ranges()->begin() ; 
    }

    do {

      for ( int ii = 0 ; ii != forvec->size() ; ii++ ) {
        if ( (*old_forvec)[ii] != (*forvec)[ii] ) {
          if ( (*forvec)[ii] == 0 ) {
            rng_maps[ii] = sub_tensops_[ii]->all_ranges()->begin() ;
          } else {
            rng_maps[ii]++;
          }
        }
      }

      old_forvec = make_shared<vector<int>>(*forvec);
      shared_ptr< vector <shared_ptr<Range_Block_Info >>> split_block = make_shared< vector <shared_ptr<Range_Block_Info >>>(num_tensors_);
      vector<shared_ptr<Range_Block_Info>>::iterator split_block_iter = split_block->begin();
      for (int jj = 0 ; jj != num_tensors_ ; jj++) 
        *split_block_iter++ = rng_maps[jj]->second;
       
      //TODO Must obtain from constraint functions 
      shared_ptr<Split_Range_Block_Info> srbi;
      {
      SRBI_Helper helper(split_block);
      srbi = make_shared<Split_Range_Block_Info>( helper, range_prime_map );
      }
      all_ranges.emplace( *(srbi->orig_block()), srbi ) ;
 
    } while( fvec_cycle_skipper( forvec, maxs, mins ) );
  

  } else { 

     //TODO  add constructor so can just have single tensor, need not be vector
    for ( auto elem : *(sub_tensops_[0]->all_ranges()) )  {
      SRBI_Helper helper( make_shared<vector<shared_ptr<Range_Block_Info >>>(1, elem.second));
      shared_ptr<Split_Range_Block_Info> srbi = make_shared<Split_Range_Block_Info>(helper, range_prime_map);
      all_ranges.emplace( *(srbi->orig_block()), srbi ) ;
    }
  }
   
  return make_shared< const map < const vector<string> , shared_ptr<Split_Range_Block_Info > > >( all_ranges );
  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::transform( const vector<char>& transformations,  const vector<int>& op_order, shared_ptr<map<char, long unsigned int>> range_prime_map  ) {  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// transformations pertain to the operators following reordering, e.g., if the origin multitensop is ABC, and the transformed version is B^{+} C A then
// op_order = { 1, 2, 0 }
// transformations = { 'H','0','0' }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::MultiTensOp<DataType>::transform" << endl;

  int num_idxs = idxs()->size();
  
  //swap round sub_tensops_ vector so can use original routines.
  vector<shared_ptr<TensOp_Base>> sub_tensops_trans = vector<shared_ptr<TensOp_Base>>(  sub_tensops_.size() );
  {
     vector<shared_ptr<TensOp_Base>>::iterator stt_it = sub_tensops_trans.begin();    
     for ( vector<int>::const_iterator  oo_it = op_order.begin(); oo_it != op_order.end(); oo_it++, ++stt_it ) 
       *stt_it = sub_tensops_[*oo_it];
  }
 
  //getting transformed idxs and aops, and performing transformations on all_ranges for sub ops  
  vector<string> trans_idxs(idxs()->size());
  vector<bool> trans_aops(aops()->size());
  {
    vector<string>::iterator ti_it = trans_idxs.begin();
    vector<bool>::iterator ta_it = trans_aops.begin();
    
    vector<char>::const_iterator t_it = transformations.begin(); 
    for ( auto stt_it = sub_tensops_trans.begin(); stt_it != sub_tensops_trans.end(); ++stt_it, t_it++  ) {
      if ( (*stt_it)->all_ranges_trans()->find(*t_it) == (*stt_it)->all_ranges_trans()->end() )
        (*stt_it)->transform(*t_it, range_prime_map);
     
      vector<string> sub_trans_idxs = vector<string>(*((*stt_it)->idxs()));
      transform_tens_vec( *t_it, sub_trans_idxs );
      copy( sub_trans_idxs.begin(), sub_trans_idxs.end(), ti_it );
      ti_it += sub_trans_idxs.size();
     
      vector<bool> sub_trans_aops = vector<bool>(*((*stt_it)->aops()));
      transform_tens_vec( *t_it, sub_trans_aops );
      copy( sub_trans_aops.begin(), sub_trans_aops.end(), ta_it );
      ta_it += sub_trans_aops.size();
    }
  }

  vector< map< const  vector<string>, shared_ptr< Range_Block_Info >>::const_iterator> rng_maps(num_tensors_);  
  vector<int> posvec( num_tensors_, 0 ); 

  map< const vector<string>, shared_ptr<Split_Range_Block_Info> > split_ranges_tmp;
  if ( num_tensors_ > 1 ) { 
    shared_ptr<vector<int>> forvec = make_shared<vector<int>>(num_tensors_, 0 ); 
    shared_ptr<vector<int>> old_forvec = make_shared<vector<int>>(*forvec);

    shared_ptr<vector<int>> mins   = make_shared<vector<int>>(num_tensors_, 0 );  
    shared_ptr<vector<int>> maxs   = make_shared<vector<int>>(num_tensors_, 0 );  

    for ( int ii = 0 ; ii != sub_tensops_trans.size() ; ii++ ){ 
      maxs->at(ii) = sub_tensops_trans[ii]->all_ranges_trans()->at(transformations[ii])->size()-1;
      rng_maps[ii] = sub_tensops_trans[ii]->all_ranges_trans()->at(transformations[ii])->begin() ; 
    }

    do {

      for ( int ii = 0 ; ii != forvec->size() ; ii++ ) {
        if ( (*old_forvec)[ii] != (*forvec)[ii] ) {
          if ( (*forvec)[ii] == 0 ) {
            rng_maps[ii] = sub_tensops_trans[ii]->all_ranges_trans()->at(transformations[ii])->begin() ;
          } else {
            rng_maps[ii]++;
          }
        }
      }

      old_forvec = make_shared<vector<int>>(*forvec);
      shared_ptr< vector <shared_ptr<Range_Block_Info >>> split_block = make_shared< vector <shared_ptr<Range_Block_Info >>>(num_tensors_);
      vector<shared_ptr<Range_Block_Info>>::iterator split_block_iter = split_block->begin();
      for (int jj = 0 ; jj != num_tensors_ ; jj++) 
        *split_block_iter++ = rng_maps[jj]->second;
       
      //TODO Must obtain from constraint functions 
      shared_ptr<Split_Range_Block_Info> srbi;
      {
      SRBI_Helper helper(split_block);
      srbi = make_shared<Split_Range_Block_Info>( helper, range_prime_map );
      }
      split_ranges_tmp.emplace( *(srbi->orig_block()), srbi ) ;
 
    } while( fvec_cycle_skipper( forvec, maxs, mins ) );
  

  } else { 

  //TODO  add constructor so can just have single tensor, need not be vector
    if ( sub_tensops_.front()->all_ranges_trans()->find(transformations.front()) == sub_tensops_.front()->all_ranges_trans()->end() )
      sub_tensops_.front()->transform(transformations.front(), range_prime_map);

    for ( auto elem : *(sub_tensops_.front()->all_ranges_trans()->at(transformations.front())) )  {
      SRBI_Helper helper( make_shared<vector<shared_ptr<Range_Block_Info >>>(1, elem.second));
      shared_ptr<Split_Range_Block_Info> srbi = make_shared<Split_Range_Block_Info>(helper, range_prime_map);
      split_ranges_tmp.emplace( *(srbi->orig_block()), srbi ) ;
    }
  }
  
  Op_dense_->split_ranges_trans()->emplace( make_pair( transformations, op_order ),  make_shared< const map < const vector<string> , shared_ptr<Split_Range_Block_Info > > >( split_ranges_tmp ) );
  
  return;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//get contractions and range, looks totally inefficient, but generating everything and then checking in sequence
//seems to work out faster; it's either check everything lots of times, or regenerate it lots of times. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void MultiTensOp::MultiTensOp<DataType>::get_ctrs_tens_ranges() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiTensOp::get_ctrs_tens_ranges " <<  endl;

  //puts uncontracted ranges into map 
  shared_ptr<vector<pair<int,int>>> noctrs = make_shared<vector<pair<int,int>>>(0);
 
  //silly, should just test {act,act,....} against contraints, and check if act in each ranges... 
  for (auto rng_it = Op_dense_->split_ranges()->begin(); rng_it != Op_dense_->split_ranges()->end(); rng_it++) {

    bool check = true;
    for (int xx = 0; xx !=rng_it->first.size() ; xx++ ) {
      if (rng_it->first[xx][0] != 'a') {
        check=false;
        break;
      }
    }

    if(!check){
      continue;

    } else {
      enter_cmtps_into_map(*noctrs, rng_it->second->factors(), rng_it->first );

    }
  }

  //puts_contractions, with specified ranges, into the map
  for ( int nctrs = 1 ; nctrs != (Op_dense_->num_idxs()/2)+1 ; nctrs++ ){
    shared_ptr<vector<shared_ptr<vector<pair<int,int>>>>> ctr_lists = get_unique_pairs( Op_dense_->plus_ops(), Op_dense_->kill_ops(), nctrs );

    for (shared_ptr<vector<pair<int,int>>> ctr_vec : *ctr_lists) {
      for (auto rng_it = Op_dense_->split_ranges()->begin(); rng_it != Op_dense_->split_ranges()->end(); rng_it++) {
  
        bool valid = true;
        //checks ranges for contractions match
        for (int ii = 0 ; ii != ctr_vec->size(); ii++){
          if ( rng_it->first[(*ctr_vec)[ii].first] != rng_it->first[(*ctr_vec)[ii].second]){
            valid = false;
            break;
          }
        }
        if (valid)
          enter_cmtps_into_map(*ctr_vec, rng_it->second->factors(), rng_it->first );
        
      }
    }
  }
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
    shared_ptr<vector<string>>  TS_id_ranges = make_shared<vector<string>>(id_ranges.begin()+Op_dense_->cmlsizevec(ii), id_ranges.begin()+cmlsizevec(ii)+sub_tensops_[ii]->num_idxs());
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
