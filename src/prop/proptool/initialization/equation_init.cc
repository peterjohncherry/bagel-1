#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/initialization/equation_init.h>
using namespace std;
using namespace bagel;
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_Value<DataType>::initialize_expressions() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_Value<DataType>::initialize_expressions()" << endl;

  cout << "range_map_ = " << endl;
  for (auto& elem : *range_map_ ) { 
    cout << "{ " <<  elem.first  << ", [ ";  cout.flush(); 
    for ( int xx : *elem.second ) 
    cout << xx << " "; cout.flush();
    cout << "] } "<< endl;
  }
  cout << endl;
  
  //TODO The looping through the terms should be on the inside, and the summation on the outside,
  //     but switching these loops round is headache, and I doubt it is significant speed wise  
  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){
    cout << "master_expression_->term_list_->at("<<ii<<").second->name_ "; cout.flush(); 
    cout << master_expression_->term_list_->at(ii).second->name_ << endl;  

    //get factors and initialize range map
    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;

    //build stuff for forvec loop in advance to simplify initialization
    shared_ptr<vector<int>> fvec = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> mins = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> maxs = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
 
    //for summations.
    shared_ptr<vector<int>> fvec_summer = make_shared<vector<int>>(*fvec);
    shared_ptr<vector<int>> mins_summer = make_shared<vector<int>>(*mins);
    shared_ptr<vector<int>> maxs_summer = make_shared<vector<int>>(*maxs) ;
 
    //build stuff for forvec loop in advance to simplify initialization
    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map; latter can have more entries
    map<string, shared_ptr<vector<int>>> term_range_map;
    vector<int>::iterator maxs_it = maxs->begin();
    for ( auto elem : *term_idx_val_map ){
      shared_ptr<vector<int>> range = range_map_->at(term_idrange_map->at(elem.first).second);
      term_range_map.emplace( elem.first, range );
      *maxs_it++ = (( range->size()-1 > 0 ) ? range->size()-1 : 0 );
    }

    vector<bool> summed_indexes(term_idrange_map->size(), false);
    {
    int qq = 0 ;
    for ( auto& elem : *term_idrange_map ) {
      
      if ( !elem.second.first ){ 
        maxs_summer->at(qq) = maxs->at(qq);
        summed_indexes[qq] = true;
      }
      qq++;
    }
    }

    do {
      // setting mins and maxs the same will prevent these indexes from being iterated in the fvec_cycler,
      // the resulting braket_list will be contain terms when all indexes but these "frozen" ones 
      // have been summed over. 
      for ( int rr = 0 ; rr != fvec->size() ;rr++ ) {
        if ( summed_indexes[rr]) {
          mins->at(rr) = fvec_summer->at(rr);
          maxs->at(rr) = fvec_summer->at(rr);
        }
      }
      copy_n(mins->begin(), mins->size(), fvec->begin());
 
      string state_ids_name = " summed:[";
      for ( auto& elem : *term_idrange_map ) 
         if (elem.second.first )
           state_ids_name+= elem.first + "," ; 
      state_ids_name += "]";

      if (state_ids_name.back() == ',' ) {
        state_ids_name.back() = ']';
      } else {
        state_ids_name += "]";
      }
  
      state_ids_name += " fixed:[";
      for ( auto& elem : *term_idrange_map ) 
         if (!elem.second.first && elem.first != "none" )
           state_ids_name+= elem.first + "," ; 
      
      if (state_ids_name.back() == ',' ) {
        state_ids_name.back() = ']';
      } else {
        state_ids_name += "]";
      }

      shared_ptr<vector<pair<DataType, string>>> expression_term_list; 
      string expression_name = master_expression_->name_ + state_ids_name;
      auto etm_loc = expression_term_map_->find( expression_name );
      if ( etm_loc == expression_term_map_->end() ) {
        expression_term_list = make_shared<vector<pair<DataType, string>>>();
        expression_term_map_->emplace( expression_name, expression_term_list );

      } else { 
        expression_term_list = etm_loc->second;
      } 

      vector<BraKet<DataType>> braket_list;;
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
              op_state_ids->at(jj) = vector<int>(bk_info.op_list_->at(jj).state_dep_);
              bk_info.op_list_->at(jj).get_op_idxs( op_state_ids->at(jj) );
            }
          }
       
          braket_list.push_back(BraKet<DataType>( bk_op_list, factor_map_->at(*bk_factors_it++), bk_info.bra_index(), bk_info.ket_index(), op_state_ids, term_init->type_ ));
        } 

      } while( fvec_cycle_skipper( fvec, maxs, mins ) );
  
      string term_name;
      int need_new_line = 0; 
      for ( BraKet<DataType>& bk : braket_list )   
        term_name  += " " + bk.bk_name() + " +" ;
      
      term_name.back()= ' ';
      term_name += state_ids_name;
      cout << endl;
      cout << master_expression_->term_list_->at(ii).second->alg_name_ << endl;
      cout << "term_name = " << term_name << endl;
      term_braket_map_->emplace( term_name, make_shared<vector<BraKet<DataType>>>(braket_list) ); 
      expression_term_list->push_back(make_pair(1.0 , term_name));

    }while( fvec_cycle_skipper( fvec_summer, maxs_summer, mins_summer ) );

  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_LinearRM<DataType>::initialize_expressions() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_LinearRM<DataType>::initialize_expression()" << endl;
  
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Equation_Init_Value<double>;
template class Equation_Init_LinearRM<double> ; 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
