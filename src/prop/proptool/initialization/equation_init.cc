#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/initialization/equation_init.h>
using namespace std;
using namespace bagel;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<MultiOp_Info>
Equation_Init_Base::get_operator_info( std::vector<std::string>& op_list, std::vector<char>& op_trans_list,
                                       std::shared_ptr<std::vector<std::vector<int>>> op_state_ids ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_Base::get_operator_info" << endl;

  print_vector( op_trans_list , " op_trans_list" ); cout << endl;

  //get state name first
  string multiop_state_name = "";
  string multiop_full_name = "";
  string multiop_name = "";
  shared_ptr<vector<shared_ptr<vector<int>>>> state_id_list =  make_shared<vector<shared_ptr<vector<int>>>>();

  shared_ptr<vector<shared_ptr<Op_Info>>> multiop_info_list = make_shared<vector<shared_ptr<Op_Info>>>( op_list.size() );
  vector<shared_ptr<Op_Info>>::iterator mil_it = multiop_info_list->begin();

  string op_name;
  string op_state_name;
  string op_full_name;
  for ( int ii = 0 ; ii != op_list.size(); ii++, mil_it++ ) {
    op_full_name = "";
    op_full_name += op_list[ii] ;

    op_name = op_full_name;

    if (op_state_ids->at(ii).size() > 0 ) {
      op_state_name +=  "_{"; 
      for( int jj = 0; jj != (*op_state_ids)[ii].size(); jj++ ) {
        op_state_name += to_string((*op_state_ids)[ii][jj]); 
      }
      op_state_name += "}"; 
    }

    op_state_name = op_full_name;

   if (op_trans_list.size() > 0 ) {
      op_full_name +=  "^{"; 
      op_full_name += op_trans_list[ii]; 
      op_full_name += "}"; 
    }
    print_vector( op_trans_list , "op_trans_list"  ) ; cout << endl;
    *mil_it = make_shared<Op_Info>( op_name, op_state_name, op_full_name,  make_shared<vector<int>> (op_state_ids->at(ii)), op_trans_list[ii] ); 

    multiop_name += op_name;
    multiop_state_name += op_state_name;
    multiop_full_name += op_full_name;
  }
  

  std::vector<int> op_order(op_list.size());
  iota( op_order.begin(), op_order.end(), 0);
  sort( op_order.begin(), op_order.end(), [ &op_list ] ( int i1, int i2) { return (bool)( op_list[i1] < op_list[i2]); });  

  return make_shared<MultiOp_Info>(multiop_name, multiop_state_name, multiop_full_name, multiop_info_list, op_order );  
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_Value<DataType>::initialize_expressions() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_Value<DataType>::initialize_expressions()" << endl;

  pair<double,double> bk_factor_dummy  = make_pair( 1.0, 1.0);
  //TODO The looping through the terms should be on the inside, and the summation on the outside,
  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){

    //get factors and initialize range map
    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;

    //build stuff for forvec loop in advance to simplify initialization
    vector<int> fvec(term_idx_val_map->size(), 0) ;
    vector<int> mins(term_idx_val_map->size(), 0) ;
    vector<int> maxs(term_idx_val_map->size(), 0) ;
 
    //for summations.
    vector<int> fvec_summer = fvec;
    vector<int> mins_summer = mins;
    vector<int> maxs_summer = maxs;
 
    //build stuff for forvec loop in advance to simplify initialization
    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map; latter can have more entries
    map<string, shared_ptr<vector<int>>> term_range_map;
    vector<string> idx_ordered_names(term_idrange_map->size()); 
    vector<int>::iterator maxs_it = maxs.begin();
    vector<string>::iterator ion_it = idx_ordered_names.begin();
    for ( auto elem : *term_idx_val_map ){
      shared_ptr<vector<int>> range = range_map_->at(term_idrange_map->at(elem.first).second);
      term_range_map.emplace( elem.first, range );
      *maxs_it++ = (( range->size()-1 > 0 ) ? range->size()-1 : 0 );
      *ion_it++ = elem.first;
    }

    vector<bool> summed_indexes(term_idrange_map->size(), false);
    {
    int qq = 0 ;
    for ( auto& elem : *term_idrange_map ) {
      if ( term_idx_val_map->find(elem.first) != term_idx_val_map->end() ) { // TODO neaten this up 
        if ( !elem.second.first ){ 
          maxs_summer[qq] = maxs[qq];
          summed_indexes[qq] = true;
        }
        qq++;
      }
    }
    }

    do {
      // setting mins and maxs the same will prevent these indexes from being iterated in the fvec_cycler,
      // the resulting braket_list will be contain terms when all indexes but these "frozen" ones 
      // have been summed over. 
      for ( int rr = 0 ; rr != fvec.size() ;rr++ ) {
        if ( summed_indexes[rr] ) {
          mins[rr] = fvec_summer[rr];
          maxs[rr] = fvec_summer[rr];
        }
      }
      copy_n(mins.begin(), mins.size(), fvec.begin());

      shared_ptr<map<string, int>> summed_map = make_shared<map<string,int>>();
      string state_ids_name = "summed : [ " ;
      for ( auto& elem : *term_idrange_map )
        if (elem.second.first )
          state_ids_name+=  " " + elem.first + ","; 
        
      if (state_ids_name.back() == ',' ) 
        state_ids_name.back() = ' ';

      state_ids_name += "]";
  
      shared_ptr<vector<pair<DataType, string>>> expression_term_list; 
      string expression_name = master_expression_->name_ + state_ids_name;
      cout << "expression_name = " << expression_name << endl;
      auto etm_loc = expression_term_map_->find( expression_name );
      if ( etm_loc == expression_term_map_->end() ) {
        expression_term_list = make_shared<vector<pair<DataType, string>>>();
        expression_term_map_->emplace( expression_name, expression_term_list );

      } else { 
        
        expression_term_list = etm_loc->second;
      } 

      vector<shared_ptr<BraKet_Base>> braket_list;
      do {
        int kk = 0 ;
        vector<int>::iterator fvec_it = fvec.begin();
        for ( auto tiv_it = term_idx_val_map->begin() ; tiv_it != term_idx_val_map->end(); tiv_it++ ) {
          if (tiv_it->first == "none" ) 
            continue ;
          tiv_it->second = term_range_map.at(tiv_it->first)->at(*fvec_it);
          fvec_it++;
        }
       
        vector<string>::iterator bk_factors_it = term_init->braket_factors_->begin();
        for ( BraKet_Init bk_info : *(term_init->braket_list_) ) {
          vector<string> bk_op_list(bk_info.op_list_->size());
          vector<char> bk_op_trans_list(bk_info.op_list_->size());
          auto op_state_ids = make_shared<vector<vector<int>>>(bk_info.op_list_->size());
          for ( int jj = 0 ; jj != bk_info.op_list_->size() ; jj++ ){
            bk_op_list[jj] = bk_info.op_list_->at(jj).name_;
            bk_op_trans_list[jj] = bk_info.op_list_->at(jj).trans_;
            if ( bk_info.op_list_->at(jj).state_dep_ > 0 ){
              op_state_ids->at(jj) = vector<int>(bk_info.op_list_->at(jj).state_dep_);
              bk_info.op_list_->at(jj).get_op_idxs( op_state_ids->at(jj) );
            }
          }

          // factor = factor_map_->at(*bk_factors_ita+)
          cout << "term_init->type_ = " << term_init->type_ << endl;
          if (term_init->type_ == "orb" ){
            shared_ptr<BraKet_OrbExcDeriv<DataType>> new_bk = make_shared< BraKet_OrbExcDeriv<DataType>>( get_operator_info( bk_op_list, bk_op_trans_list, op_state_ids ), bk_factor_dummy,
                                                                                              bk_info.bra_index(), bk_info.ket_index(), term_init->type_);
            new_bk->target_op_ = term_init->proj_op_name_; 
            cout << "new_bk.target_op_  = " << new_bk->target_op_  << endl;
            cout << "term_init->proj_op_name_ = " <<  term_init->proj_op_name_ << endl; 
            new_bk->orb_exc_deriv_ = true;
            braket_list.push_back( new_bk ) ;
          } else {
            braket_list.push_back(make_shared<BraKet_Full<DataType>>( get_operator_info( bk_op_list, bk_op_trans_list, op_state_ids ), bk_factor_dummy,
                                                    bk_info.bra_index(), bk_info.ket_index(), term_init->type_));
          }

        } 

      } while( fvec_cycle_skipper( fvec, maxs, mins ) );
       
      vector<pair<string,int>> fixed_id_vals; 
      state_ids_name += " fixed:[";
      for ( int rr = 0 ; rr != fvec.size() ;rr++ ) {
        if ( summed_indexes[rr] && idx_ordered_names[rr] != "none" ){ 
	   state_ids_name += "( " + idx_ordered_names[rr] + " : " + to_string(term_idx_val_map->at(idx_ordered_names[rr])) +" ),"; 
           fixed_id_vals.push_back(make_pair(idx_ordered_names[rr], term_idx_val_map->at(idx_ordered_names[rr]) ));
        }
      }

      if (state_ids_name.back() == ',' )
        state_ids_name.back() = ' ';
      state_ids_name += "]";
 
      string term_name;
      int need_new_line = 0; 
      for ( shared_ptr<BraKet_Base>& bk : braket_list )   
        term_name  += " " + bk->bk_name() + " +" ;
      
      term_name.back()= ' ';
      term_name += state_ids_name;
      cout << endl;
      cout << "term_name = " << term_name << endl;
      term_braket_map_->emplace( term_name, make_shared<vector<shared_ptr<BraKet_Base>>>(braket_list) ); 
      expression_term_list->push_back( make_pair( (DataType)1.0 , term_name) );
       
      expression_term_map_->emplace( term_name , make_shared<vector<pair<DataType, string>>>( 1,  make_pair( (DataType)1.0 , term_name) ) );

      string term_alg_name = master_expression_->term_list_->at(ii).second->alg_name_;;
       
      expression_term_map_by_states_->emplace( make_pair(term_alg_name, fixed_id_vals), make_shared<vector<pair<DataType, string>>>( 1,  make_pair( (DataType)1.0 , term_name) ) );

    } while( fvec_cycle_skipper( fvec_summer, maxs_summer, mins_summer ) );

  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Equation_Init_Value<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
