#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/initialization/equation_init.h>
using namespace std;
using namespace bagel;
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_Value<DataType>::initialize_expressions() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_EQUATION_INIT_VALUE 
cout << "Equation_Init_Value<DataType>::initialize_expressions()" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////

  pair<double,double> bk_factor_dummy  = make_pair( 1.0, 0.0);
  //TODO The looping through the terms should be on the inside, and the summation on the outside,
  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){

    //get factors and initialize range map
    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;

    //build stuff for forvec loop in advance to simplify initialization
    vector<int> fvec_fixed_ids(term_idx_val_map->size(), 0) ;
    vector<int> mins_fixed_ids(term_idx_val_map->size(), 0) ;
    vector<int> maxs_fixed_ids(term_idx_val_map->size()) ;
 
    //build stuff for forvec loop in advance to simplify initialization
    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map; latter can have more entries
    map<string, shared_ptr<vector<int>>> term_range_map;
    vector<string> idx_ordered_names(term_idrange_map->size()); 
    vector<int>::iterator mfi_it = maxs_fixed_ids.begin();
    vector<string>::iterator ion_it = idx_ordered_names.begin();

    for ( auto elem : *term_idx_val_map ){
      shared_ptr<vector<int>> range = range_map_->at(term_idrange_map->at(elem.first).second);
      term_range_map.emplace( elem.first, range );
      *mfi_it++ = ( range->size()-1 > 0 ) ? range->size()-1 : 0 ;
      *ion_it++ = elem.first;
    }

    vector<bool> fixed_indexes( term_idx_val_map->size(), false ) ;
    {
    int qq = 0 ;
    for ( auto& elem : *term_idrange_map ) {
      if ( term_idx_val_map->find(elem.first) != term_idx_val_map->end() ) { // TODO neaten this up 
        if ( !elem.second.first ){ 
          fixed_indexes[qq] = true;
        } else {
          fixed_indexes[qq] = false;
          fvec_fixed_ids[qq] = maxs_fixed_ids[qq]; 
          mins_fixed_ids[qq] = maxs_fixed_ids[qq];
        }
        qq++;
      }
    }
    }
    string summed_ids_name = "summed_ids : [ " ;
    for ( auto& elem : *term_idrange_map )
      if (elem.second.first )
        summed_ids_name+=  " " + elem.first + ","; 
    summed_ids_name.back() = ' ';
    summed_ids_name += ']';

    do {
      vector<int> fvec_summed_ids( mins_fixed_ids.size() , 0 );
      vector<int> mins_summed_ids( mins_fixed_ids.size() , 0 );
      vector<int> maxs_summed_ids = maxs_fixed_ids;

      vector<pair<string,int>> fixed_id_vals; 
      for ( int rr = 0 ; rr != fixed_indexes.size() ;rr++ ){
        if ( fixed_indexes[rr] ){
          fixed_id_vals.push_back(make_pair(idx_ordered_names[rr], term_idx_val_map->at(idx_ordered_names[rr]) ));
          fvec_summed_ids[rr] = fvec_fixed_ids[rr];
          mins_summed_ids[rr] = fvec_fixed_ids[rr];
          maxs_summed_ids[rr] = fvec_fixed_ids[rr];
        }
      }

      //for summations.
      vector<shared_ptr<BraKet_Base>> braket_list(0);
      do {
        vector<int>::iterator fsi_it = fvec_summed_ids.begin();
        for ( auto tiv_it = term_idx_val_map->begin() ; tiv_it != term_idx_val_map->end(); tiv_it++,  fsi_it++ ) {
          if (tiv_it->first == "none" ) 
            continue ;
          tiv_it->second = term_range_map.at(tiv_it->first)->at(*fsi_it);
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

          shared_ptr<MultiOp_Info> multiop_info =  make_shared<MultiOp_Info>( bk_op_list, bk_op_trans_list, op_state_ids ); 
          if ( multiop_info->canonical_order_ )
            multiop_info->op_info_canonical_ = dynamic_pointer_cast<MultiOp_Info>(multiop_info->shared_from_this());

          if (term_init->type_ == "orb" ){
            shared_ptr<BraKet_OrbExcDeriv<DataType>> new_bk = make_shared< BraKet_OrbExcDeriv<DataType>>( multiop_info, bk_factor_dummy, bk_info.bra_index(), bk_info.ket_index(), term_init->type_);
            new_bk->target_op_ = term_init->proj_op_name_; 
            new_bk->orb_exc_deriv_ = true;
            braket_list.push_back( new_bk ) ;
          } else {
            braket_list.push_back(make_shared<BraKet_Full<DataType>>( multiop_info, bk_factor_dummy, bk_info.bra_index(), bk_info.ket_index(), term_init->type_));
          }

        }
      } while( fvec_cycle_skipper( fvec_summed_ids, maxs_summed_ids, mins_summed_ids ) );

      string state_ids_name = summed_ids_name;
      state_ids_name += " fixed_ids : [ ";  
      for ( int rr = 0 ; rr != fvec_fixed_ids.size() ;rr++ )
        if ( fixed_indexes[rr]  )
          state_ids_name += "( " + idx_ordered_names[rr] + " : " + to_string(term_idx_val_map->at(idx_ordered_names[rr])) +" ),"; 
      state_ids_name.back() = ']';

      string term_name =  state_ids_name;
      for ( shared_ptr<BraKet_Base>& bk : braket_list )   
        term_name  += " " + bk->bk_name() + " +" ;
      term_name.back() = ' ';

      string expression_name = master_expression_->name_ +  state_ids_name;
      shared_ptr<vector<pair<DataType, string>>> expression_term_list; 
      auto etm_loc = expression_term_map_->find( expression_name );
      if ( etm_loc == expression_term_map_->end() ) {
        expression_term_list = make_shared<vector<pair<DataType, string>>>(0);
        expression_term_map_->emplace( expression_name, expression_term_list );
      } else {
        expression_term_list = etm_loc->second;
      } 

      expression_term_list->push_back( make_pair( (DataType)1.0 , term_name) );

      term_braket_map_->emplace( term_name, make_shared<vector<shared_ptr<BraKet_Base>>>(braket_list) );

      string term_alg_name = master_expression_->term_list_->at(ii).second->alg_name_;

      // expression_term_map_->emplace( term_name , make_shared<vector<pair<DataType, string>>>( 1,  make_pair( (DataType)1.0 , term_name) ) );
      expression_term_map_by_states_->emplace( make_pair(term_alg_name, fixed_id_vals), make_shared<vector<pair<DataType, string>>>( 1,  make_pair( (DataType)1.0 , term_name) ) );
    } while( fvec_cycle_skipper( fvec_fixed_ids, maxs_fixed_ids, mins_fixed_ids ) );

  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Equation_Init_Value<double>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
