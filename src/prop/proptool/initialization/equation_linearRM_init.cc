#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/initialization/equation_init.h>
using namespace std;
using namespace bagel;
/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType> 
void Equation_Init_LinearRM<DataType>::initialize_expressions() {
/////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Init_LinearRM<DataType>::initialize_expressions()" << endl;

  for ( int ii = 0 ; ii != master_expression_->term_list_->size(); ii++  ){
    cout << " ii = " << ii << endl;
    int counter = 0;
    DataType term_factor = factor_map_->at(master_expression_->term_list_->at(ii).first);
    shared_ptr< Term_Init > term_init = master_expression_->term_list_->at(ii).second;
    expression_term_map_->emplace( term_init->name_ , make_shared<vector<pair<DataType, string>>>( 1,  make_pair( (DataType)1.0 , term_init->name_) ) );
    shared_ptr<map< string, pair<bool, string>>> term_idrange_map = master_expression_->term_range_maps_->at(ii);
    shared_ptr<map< string, int>> term_idx_val_map = term_init->idx_val_map_;

    // build stuff for forvec loop in advance to simplify initialization
    shared_ptr<vector<int>> fvec = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> mins = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
    shared_ptr<vector<int>> maxs = make_shared<vector<int>>(term_idx_val_map->size(), 0) ;
 
    // build stuff for forvec loop in advance to simplify initialization
    // NOTE : because ordering of term_init->idx_val_map is not necessarily same as term_idrange_map; latter can have more entries
    map< string, shared_ptr<vector<int>> > term_range_map;
    vector<string> idx_ordered_names( term_idrange_map->size() ); 
    vector<int>::iterator maxs_it = maxs->begin();
    vector<string>::iterator ion_it = idx_ordered_names.begin();
    for ( auto elem : *term_idx_val_map ){
      shared_ptr<vector<int>> range = range_map_->at(term_idrange_map->at(elem.first).second);
      term_range_map.emplace( elem.first, range );
      *maxs_it++ = (( range->size()-1 > 0 ) ? range->size()-1 : 0 );
      *ion_it++ = elem.first;
    }
  
    do {

      //print_vector( *fvec, "fvec" ) ; cout.flush();
      int kk = 0 ;
      vector<int>::iterator fvec_it = fvec->begin();
      for ( auto tiv_it = term_idx_val_map->begin() ; tiv_it != term_idx_val_map->end(); tiv_it++ ) {
        if (tiv_it->first == "none" ) { continue ; }
        tiv_it->second = term_range_map.at(tiv_it->first)->at(*fvec_it);
        fvec_it++;
      }
      
      vector<BraKet<DataType>> braket_list;
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
      
      vector<pair<string,int>> fixed_id_vals; 

      cout << term_init->name_ << " " ; cout.flush(); cout << " [ "; cout.flush();  
      for ( int rr = 0 ; rr != fvec->size() ;rr++ ) {  
        fixed_id_vals.push_back(make_pair( idx_ordered_names[rr], term_idx_val_map->at(idx_ordered_names[rr]) ));
        cout << "{ " << fixed_id_vals.back().first << ", " << fixed_id_vals.back().second << " } "; cout.flush();
      }
      cout <<  "] " <<  endl;
      sort(fixed_id_vals.begin(), fixed_id_vals.end()); 

      term_braket_map_->emplace( term_init->name_, make_shared<vector<BraKet<DataType>>>(braket_list) ); 

      expression_term_map_by_states_->emplace( make_pair(term_init->name_, fixed_id_vals), make_shared<vector<pair<DataType, string>>>( 1, make_pair((DataType)1.0, term_init->name_)));

    } while( fvec_cycle_skipper( fvec, maxs, mins ) );

  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Equation_Init_LinearRM<double> ; 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
