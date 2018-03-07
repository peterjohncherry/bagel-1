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
void BraKet<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<TensOp_Base>>> MT_map,                
                                                            shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo> >>>> G_to_A_map,
                                                            shared_ptr<map<string, shared_ptr< GammaInfo >>> gamma_info_map,
                                                            shared_ptr<StatesInfo<DataType>> target_states,
                                                            shared_ptr<set<string>> required_blocks ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::generate_gamma_Atensor_contractions : " << name_ << endl; 
  //TODO fix this so it uses proper number of states; if statement in center should call Bra_num Ket_num appropriate Ops
  //     GammaGen should be initialized outside states loop, and wipe Gamma_Vec for each new range.
  //     Loop through dense ranges on the outside, then check sparsity on the inner when adding to gamma_info_map.
  
  Total_Op_ = MT_map->at(multiop_name_);

  vector<char> projector_names; 
  for ( auto& tens : Total_Op_->sub_tensops() ) 
    if ( tens->is_projector() ) 
      projector_names.push_back(tens->name()[0]);

  shared_ptr<vector<string>> idxs_buff  = make_shared<vector<string>>(*(Total_Op_->idxs()));
  shared_ptr<vector<bool>> aops_buff  = make_shared<vector<bool>>(*Total_Op_->aops());    

  for ( auto range_map_it = Total_Op_->split_ranges()->begin(); range_map_it !=Total_Op_->split_ranges()->end(); range_map_it++ ){
    if ( range_map_it->second->survives() && !range_map_it->second->is_sparse( op_state_ids_ ) ){  

      auto GGen = make_shared<GammaGenerator>( target_states, bra_num_, ket_num_, idxs_buff, aops_buff, gamma_info_map, G_to_A_map, factor_ );
      GGen->add_gamma( range_map_it->second );
      
      if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){
        if ( GGen->generic_reorderer( "normal order", false, false ) ) {
          if ( GGen->generic_reorderer( "alternating order", false, true ) ){  

            cout << "We need these blocks : " ; cout.flush(); cout << " Total_Op_->sub_tensops().size() = " ; cout.flush(); cout << Total_Op_->sub_tensops().size() << endl; 
            vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();
            int qq = 0 ;
           for ( auto&  tens_block : *(range_map_it->second->range_blocks()) ){ 
              cout << tens_block->orig_name()  << " " ; cout.flush(); cout << "sub_tensops[" << qq<< " ]->name() = "; cout.flush(); cout << sub_tensops[qq]->name() << endl;
              MT_map->at( sub_tensops[qq++]->name()  )->add_required_block( tens_block->orig_name() );
              required_blocks->emplace( tens_block->orig_name() );
            }
            cout << endl; 
          }
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
        for ( int qq = 0; qq !=A_map_it->second->id_orders().size(); qq++) {
          cout << "{ " ; print_vector( A_map_it->second->id_order(qq), "" ); 
          cout << " (" << A_map_it->second->factor(qq).first <<  "," <<  A_map_it->second->factor(qq).first<<  ") }" ; cout <<endl;
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
