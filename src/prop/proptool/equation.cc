#include <src/prop/proptool/equation.h>
#include <src/smith/wicktool/wickutils.h>
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
    
    do {
   
      vector<int>::iterator fvec_it = fvec->begin(); 
      for ( auto elem : *term_idx_val_map )
        elem.second = term_range_map.at(elem.first)->at(*fvec_it++);
   
      for ( BraKet_Init bk_info : *(term_init->braket_list_) ) {
        cout << "<" << bk_info.bra_index() << "|" << endl; 
        auto op_idxs_list =  bk_info.get_op_idxs_list(); 
        for ( int jj = 0 ; jj != bk_info.op_list_->size() ; jj++ ){ 
          cout << bk_info.op_list_->at(jj).name_<< "_{ "; 
          for ( int op_id : op_idxs_list->at(jj) )
            cout << op_id ; 
          cout << "}" ;
        } 
        cout << "|" << bk_info.ket_index() << ">" << endl; 
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
     cout << "ii = " << ii << endl; 
    //get factors and initialize range map    ;
    cout << "factor_map_->at( " << master_expression_->term_list_->at(ii).first<< ") = " ; cout.flush(); 
    cout <<  factor_map_->at(master_expression_->term_list_->at(ii).first) <<  endl;

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
    for ( auto elem : *term_idx_val_map ){ 
      *maxs_it = term_range_map.at(elem.first)->size();
       if ( *maxs_it != 0 )
         (*maxs_it) -=1;
       maxs_it++;
    }
    do {
   
    print_vector( *fvec , "fvec") ; cout << endl;
      int kk = 0 ;
      for ( auto elem : *term_idx_val_map ){
        cout << "elem.first = " << elem.first ; cout.flush(); cout << " kk = " << kk << "   term_range_map.at(" << elem.first << ")->at(" << fvec->at(kk)<< " );" ; cout.flush();
        if (elem.first  == "none" ) { cout << endl;  continue; } 
        cout << term_range_map.at(elem.first)->at(fvec->at(kk)) ; cout.flush();
        cout << "    elem.second = " << elem.second << endl;
        elem.second = term_range_map.at(elem.first)->at(fvec->at(kk));
        kk++;
      }
      for ( BraKet_Init bk_info : *(term_init->braket_list_) ) {
        cout << "<" << bk_info.bra_index() << "|" << endl; 
        auto op_idxs_list =  bk_info.get_op_idxs_list(); 
        for ( int jj = 0 ; jj != bk_info.op_list_->size() ; jj++ ){ 
          cout << bk_info.op_list_->at(jj).name_<< "_{ "; 
          for ( int op_id : op_idxs_list->at(jj) )
            cout << op_id ; 
          cout << "}" ;
        } 
        cout << "|" << bk_info.ket_index() << ">" << endl; 
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


