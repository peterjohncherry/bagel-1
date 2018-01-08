#include <bagel_config.h>
#ifdef COMPILE_SMITH
  #include <src/smith/wicktool/braket.h>
 //#include "braket.h"

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::Build_TotalOp(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::Build_TotalOp" << endl;
  string MT_name = "";
  for( shared_ptr<TensOp::TensOp<DataType>> Tens : *Sub_Ops )
    MT_name += Tens->name();

  Total_Op = make_shared<MultiTensOp::MultiTensOp<DataType>>( MT_name , true, *Sub_Ops );
  Total_Op->get_ctrs_tens_ranges();

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::Build_Gamma_SpinFree(shared_ptr<const vector<bool>> aops, shared_ptr<const vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::Build_Gamma_SpinFree_New" << endl;
  //TODO fix this so it uses proper number of states; if statement in center should call Bra_num Ket_num appropriate Ops
  //     GammaGen should be initialized outside states loop, and wipe Gamma_Vec for each new range.
  //     Loop through dense ranges on the outside, then check sparsity on the inner when adding to GammaMap.
  int nstates = 1;
 
  shared_ptr<vector<string>> idxs_buff  = make_shared<vector<string>>(*idxs );
  shared_ptr<vector<bool>> aops_buff  = make_shared<vector<bool>>(*aops );    
  
  for ( int Ket_num = 0 ; Ket_num != nstates; Ket_num++ ){
    for ( int Bra_num = 0 ; Bra_num != nstates; Bra_num++ ){
      for ( auto range_map_it = Total_Op->all_ranges()->begin(); range_map_it !=Total_Op->all_ranges()->end(); range_map_it++ ){
        
        print_vector( *(range_map_it->second->unique_block()), " ranges into gamma? " );  if ( !range_map_it->second->survives() ) {  cout << " ... no " << endl; } 

        if ( range_map_it->second->survives() ) {  cout << " ... yes " << endl; 
          shared_ptr<GammaGenerator>  GGen = make_shared<GammaGenerator>(TargetStates, Bra_num, Ket_num, idxs_buff, aops_buff, GammaMap, G_to_A_map, factor_); 
          GGen->add_gamma( range_map_it->second );
          GGen->norm_order();
          GGen->optimized_alt_order();
        } 
      }
    }
  }
 
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){

    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;

    if ( map_it->second->size() == 0 ) { 
       
      cout << endl << "      - No Contributions! -" <<  endl << endl;    
          
    } else {

      for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
        cout <<  A_map_it->first << "  "; cout.flush();
        for ( int qq = 0; qq !=A_map_it->second.id_orders.size(); qq++) {
          cout << "{ " ; print_vector( A_map_it->second.id_order(qq), "" ); 
          cout << " (" << A_map_it->second.factor(qq).first <<  "," <<  A_map_it->second.factor(qq).first<<  ") }" ; cout <<endl;
        }
        cout << endl;
      }
    }
  }

  return; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::Build_Gamma_WithSpin(shared_ptr<const vector<bool>> aops, shared_ptr< const vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_WithSpin" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pairs up each gamma term (defined by a spin pathway,and spin sector)  with relevant CtrMultiTensorParts .
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::Build_Tensor_Contraction_list_CMTP(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
