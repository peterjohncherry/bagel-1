#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
 //#include "braket.h"

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<MultiTensOp::MultiTensOp<DataType>>>> MT_map,                
                                                            shared_ptr<map<string, shared_ptr< map<string, AContribInfo >>>> G_to_A_map,
                                                            shared_ptr<map<string, shared_ptr< GammaInfo >>> gamma_info_map,
                                                            shared_ptr<StatesInfo<DataType>> target_states ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::generate_gamma_Atensor_contractions : " << name_ << endl; 
  //TODO fix this so it uses proper number of states; if statement in center should call Bra_num Ket_num appropriate Ops
  //     GammaGen should be initialized outside states loop, and wipe Gamma_Vec for each new range.
  //     Loop through dense ranges on the outside, then check sparsity on the inner when adding to gamma_info_map.
  
  Total_Op_ = MT_map->at(multiop_name_);
  shared_ptr<vector<string>> idxs_buff  = make_shared<vector<string>>(*(Total_Op_->idxs()));
  shared_ptr<vector<bool>> aops_buff  = make_shared<vector<bool>>(*Total_Op_->aops());    

  for ( auto range_map_it = Total_Op_->split_ranges()->begin(); range_map_it !=Total_Op_->split_ranges()->end(); range_map_it++ ){
    
    print_vector( *(range_map_it->second->unique_block()), " ranges into gamma? " ); 

    if ( !range_map_it->second->survives() ) {
      cout << " ... no ; needs to be contracted with something else. " << endl;  

    } else if ( range_map_it->second->is_sparse( op_state_ids_ ) ){  
      cout << " ... no ; this block is sparse" << endl;  

    } else { 
      cout << " ... yes " << endl; 
      
      shared_ptr<GammaGenerator>  GGen = make_shared<GammaGenerator>( target_states, bra_num_, ket_num_, idxs_buff, aops_buff, gamma_info_map, G_to_A_map, factor_ );
      GGen->add_gamma( range_map_it->second );
      GGen->norm_order();
      bool does_this_block_contribute = GGen->optimized_alt_order();
      if ( does_this_block_contribute ) {
        cout << "We need these blocks : " ; cout.flush();
        auto tvec_it = Total_Op_->orig_tensors_.begin();
        for ( auto&  tens_block : *(range_map_it->second->range_blocks()) ){ 
          cout << tens_block->orig_name()  << " " ; cout.flush(); 
          cout << "(*tvec_it)->name() = " << (*tvec_it)->name() << endl;
          (*tvec_it)->add_required_block( tens_block->orig_name() );
          tvec_it++; 
        }
        cout << endl; 
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<double>;
//template class BraKet<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
