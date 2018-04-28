#include <src/prop/proptool/algebraic_manipulator/equation_linearRM.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
// Terms are state specific parts of the expression
//////////////////////////////////////////////////////////////////////////
using namespace std;
template<typename DataType>
void Equation_LinearRM<DataType>::generate_state_specific_terms() {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_LinearRM<DataType>::generate_state_specific_terms() " << endl;  

  for ( auto term_info : *term_braket_map_state_spec_ ){
   
    vector<pair<string,int>> fixed_idxs = term_info.first.second;
    sort( fixed_idxs.begin(), fixed_idxs.end() );

    cout << term_info.first.first << ": [ ";
    for ( auto  fi_it = fixed_idxs.begin() ; fi_it != fixed_idxs.end(); fi_it++ ) 
      cout << "(" << fi_it->first << "," << fi_it->second << ") " ;
    cout << "]" << endl;

    for ( auto  fi_it = fixed_idxs.begin() ; fi_it != fixed_idxs.end(); fi_it++ ) 
      if( fi_it->first == "none") {  
        fixed_idxs.erase( fi_it ); 
        break;
      }
    
    cout << "fixed_indexes : " <<  endl;
    for ( auto  fi_it = fixed_idxs.begin() ; fi_it != fixed_idxs.end(); fi_it++ ) 
      cout << "(" << fi_it->first << "," << fi_it->second << ") " ;
    cout << endl;

    auto new_key = make_pair( term_info.first.first, fixed_idxs );  
    if ( term_map_->find( new_key ) == term_map_->end() ) 
      add_term( new_key, term_info.second );

  }
  cout << " LEAVING Equation_LinearRM<DataType>::generate_state_specific_terms() " << endl;  
  return;

}

//////////////////////////////////////////////////////////////////////////
// Terms are state specific parts of the expression
//////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
Equation_LinearRM<DataType>::add_term( pair<string, vector<pair<string,int>>>&  new_key,
                                       shared_ptr<vector<std::shared_ptr<BraKet_Base>>> expr_bk_list ) {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_LinearRM<DataType>::add_term" << endl;

  string expression_type = this->add_expression_info( expr_bk_list );

  if ( expression_type == "orbital_excitation_derivative"  ) {
    term_map_->emplace( new_key, make_shared<Expression_Orb_Exc_Deriv<DataType>>( expr_bk_list, states_info_, MT_map_, CTP_map_, ACompute_map_, gamma_info_map_, expression_type ));
  } else  if ( expression_type == "full"  ) {
    term_map_->emplace( new_key, make_shared<Expression_Full<DataType>>( expr_bk_list, states_info_, MT_map_, CTP_map_, ACompute_map_, gamma_info_map_, expression_type ));
  } else {
    assert (false);
//    throw std::logic_error( "have not implemented expression type \"" + expression_type "\" ... Aborting!!" );
  }
  return;
}
//////////////////////////////////////////////////////////////////////////
template class Equation_LinearRM<double>;
////////////////////////////////////////////////////////////////////////// 
