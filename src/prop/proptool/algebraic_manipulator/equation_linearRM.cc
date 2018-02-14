#include <src/prop/proptool/algebraic_manipulator/equation_linearRM.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////
// Terms are state specific parts of the expression
//////////////////////////////////////////////////////////////////////////
using namespace std;
template<typename DataType>
void Equation_LinearRM<DataType>::generate_state_specific_terms() {  
//////////////////////////////////////////////////////////////////////////
cout << " void Equation_LinearRM<DataType>::generate_all_terms() " << endl;  
  
  for ( auto& term_info : *term_braket_map_state_spec_ ) 
    if ( term_map_->find( term_info.first ) == term_map_->end() )
      term_map_->emplace( term_info.first, this->build_expression( term_info.second ) );

  return;

}
//////////////////////////////////////////////////////////////////////////
template class Equation_LinearRM<double>;
////////////////////////////////////////////////////////////////////////// 
