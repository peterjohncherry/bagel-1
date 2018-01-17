#include <bagel_config.h>
#ifdef COMPILE_SMITH
  #include <src/smith/wicktool/braket.h>
 //#include "braket.h"

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
BraKet<DataType>::BraKet( BraKet_Init<DataType>& bk_info, 
                          shared_ptr<map<string,shared_ptr<MultiTensOp::MultiTensOp<DataType>>>> MT_map,                
                          shared_ptr<map<string, shared_ptr< map<string, AContribInfo >>>> G_to_A_map,
                          shared_ptr<map<string, shared_ptr< GammaInfo >>> GammaMap,
                          shared_ptr<StatesInfo<DataType>> target_states ) : 
                          Total_Op_(MT_map->at(bk_info.multiop)),
                          bra_num_(bk_info.bra_num_), ket_num_(bk_info.ket_num_), factor_(bk_info.factor) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::BraKet" << endl;
  //TODO fix this so it uses proper number of states; if statement in center should call Bra_num Ket_num appropriate Ops
  //     GammaGen should be initialized outside states loop, and wipe Gamma_Vec for each new range.
  //     Loop through dense ranges on the outside, then check sparsity on the inner when adding to GammaMap.
  
  shared_ptr<vector<string>> idxs_buff  = make_shared<vector<string>>(*(Total_Op_->idxs()));
  shared_ptr<vector<bool>> aops_buff  = make_shared<vector<bool>>(*Total_Op_->aops());    

  for ( auto range_map_it = Total_Op_->all_ranges()->begin(); range_map_it !=Total_Op_->all_ranges()->end(); range_map_it++ ){
    
    print_vector( *(range_map_it->second->unique_block()), " ranges into gamma? " );  if ( !range_map_it->second->survives() ) {  cout << " ... no " << endl; } 

    if ( range_map_it->second->survives() ) {  cout << " ... yes " << endl; 
      shared_ptr<GammaGenerator>  GGen = make_shared<GammaGenerator>( target_states, bra_num_, ket_num_, idxs_buff, aops_buff, GammaMap, G_to_A_map, factor_ ); 
      GGen->add_gamma( range_map_it->second );
      GGen->norm_order();
      bool does_this_block_contribute = GGen->optimized_alt_order();

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
