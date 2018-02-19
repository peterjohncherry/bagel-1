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
    
    cout << "without none" <<  endl;
    for ( auto  fi_it = fixed_idxs.begin() ; fi_it != fixed_idxs.end(); fi_it++ ) 
      cout << "(" << fi_it->first << "," << fi_it->second << ") " ;

    auto new_key = make_pair( term_info.first.first, fixed_idxs );  
    cout <<"asdxx1"<<endl;
    if ( term_map_->find( new_key ) == term_map_->end() ) 
      term_map_->emplace( new_key, this->build_expression( term_info.second ) );
    cout <<"asdxx2"<<endl;

  }
  cout << " LEAVING Equation_LinearRM<DataType>::generate_state_specific_terms() " << endl;  
  return;

}
//////////////////////////////////////////////////////////////////////////
template class Equation_LinearRM<double>;
////////////////////////////////////////////////////////////////////////// 
