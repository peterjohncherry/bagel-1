#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/initialization/equation_init.h>
using namespace std;
using namespace bagel;
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_Value<DataType>::initialize_expression() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_Value<DataType>::initialize_expression()" << endl;
   
  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){
   
    //get factors and initialize range map
    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;

    //build stuff for forvec loop in advance to simplify initialization
    shared_ptr<vector<int>> fvec = make_shared<vector<int>>(term_idrange_map->size(), 0) ;
    shared_ptr<vector<int>> mins = make_shared<vector<int>>(term_idrange_map->size(), 0) ;
    shared_ptr<vector<int>> maxs = make_shared<vector<int>>(term_idrange_map->size(), 0) ;
 
    map<string, shared_ptr<vector<int>>> term_range_map;
    for ( auto elem : *term_idrange_map ) {
      shared_ptr<vector<int>> range = range_map_->at(elem.second.second);
      term_range_map.emplace( elem.first, range);
    }

    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map;
    vector<int>::iterator maxs_it = maxs->begin();
    for ( auto elem : *term_idx_val_map )
      *maxs_it++ = term_range_map.at(elem.first)->size()-1;

 
    vector<BraKet<DataType>> braket_list; 
    do {
      int kk = 0 ;
      vector<int>::iterator fvec_it = fvec->begin();
      for ( auto tiv_it = term_idx_val_map->begin() ; tiv_it != term_idx_val_map->end(); tiv_it++ ) {
         if (tiv_it->first == "none" ) { continue ; }
         tiv_it->second = term_range_map.at(tiv_it->first)->at(*fvec_it);
         fvec_it++;
      }

      vector<string>::iterator bk_factors_it = term_init->braket_factors_->begin();
      for ( BraKet_Init bk_info : *(term_init->braket_list_) ) {
        vector<string> bk_op_list(bk_info.op_list_->size());
        auto op_state_ids = make_shared<vector<vector<int>>>(bk_info.op_list_->size());
        for ( int jj = 0 ; jj != bk_info.op_list_->size() ; jj++ ){
          bk_op_list[jj] = bk_info.op_list_->at(jj).name_;
          if ( bk_info.op_list_->at(jj).state_dep_ > 0 ){
            op_state_ids->at(jj) =  vector<int>(bk_info.op_list_->at(jj).state_dep_);
            bk_info.op_list_->at(jj).get_op_idxs( op_state_ids->at(jj) ) ;
          }
        }

        braket_list.push_back(BraKet<DataType>( bk_op_list, factor_map_->at(*bk_factors_it++), bk_info.bra_index(), bk_info.ket_index(), op_state_ids, term_init->type_ ));
     } 
 
    } while(fvec_cycle_skipper( fvec, maxs, mins )) ;

    // loop through vectors in range map, to change what is in id val map
    //
    // on each loop, go through braket_list in term init, and get the op indexes and bra and ket indexes
    //
    // use these outputs to initialize the BraKet list, which can then be fed into Expression. 
    //
    // emplace this expression into the map.

  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_LinearRM<DataType>::initialize_expression() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_Value<DataType>::initialize_expression()" << endl;
  
  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){ 

    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;  

    //build stuff for forvec loop in advance to simplify initialization
    shared_ptr<vector<int>> fvec = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> mins = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> maxs = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
 
    map<string, shared_ptr<vector<int>>> term_range_map;
    for ( auto elem : *term_idrange_map ) {
      shared_ptr<vector<int>> range = range_map_->at(elem.second.second);   
      term_range_map.emplace( elem.first, range);
    }
    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map;
    vector<int>::iterator maxs_it = maxs->begin();
    vector<string> name_list(maxs->size()); 
    vector<string>::iterator name_list_it  = name_list.begin();
    for ( auto elem : *term_idx_val_map ){ 
      *maxs_it = term_range_map.at(elem.first)->size();
      *name_list_it++ = elem.first; 
       if ( *maxs_it != 0 )
         (*maxs_it) -=1;
       maxs_it++;
       
    }
    
    vector<BraKet<DataType>> braket_list; 
    do {
      int kk = 0 ;
      vector<int>::iterator fvec_it = fvec->begin();
      for ( auto tiv_it = term_idx_val_map->begin() ; tiv_it != term_idx_val_map->end(); tiv_it++ ) {
         if (tiv_it->first == "none" ) { continue ; }
         tiv_it->second = term_range_map.at(tiv_it->first)->at(*fvec_it);
         fvec_it++;
      }

      vector<string>::iterator bk_factors_it = term_init->braket_factors_->begin();
      for ( BraKet_Init bk_info : *(term_init->braket_list_) ) {
        vector<string> bk_op_list(bk_info.op_list_->size());
        auto op_state_ids = make_shared<vector<vector<int>>>(bk_info.op_list_->size());
        for ( int jj = 0 ; jj != bk_info.op_list_->size() ; jj++ ){
          bk_op_list[jj] = bk_info.op_list_->at(jj).name_;
          if ( bk_info.op_list_->at(jj).state_dep_ > 0 ){
            op_state_ids->at(jj) =  vector<int>(bk_info.op_list_->at(jj).state_dep_);
            bk_info.op_list_->at(jj).get_op_idxs( op_state_ids->at(jj) ) ;
          }
        }
        
        BraKet<DataType>  new_bk( bk_op_list, factor_map_->at(*bk_factors_it++), bk_info.bra_index(), bk_info.ket_index(), op_state_ids, term_init->type_ );
        cout <<  new_bk.bk_name() << endl;
        braket_list.push_back( new_bk );
      } 
 
    } while(fvec_cycle_skipper( fvec, maxs, mins )) ;

    // loop through vectors in range map, to change what is in id val map
    //
    // on each loop, go through braket_list in term init, and get the op indexes and bra and ket indexes
    //
    // use these outputs to initialize the BraKet list, which can then be fed into Expression. 
    //
    // emplace this expression into the map.

  }
  return;

}


template class Equation_Init_Value<double>;
template class Equation_Init_LinearRM<double> ; 


