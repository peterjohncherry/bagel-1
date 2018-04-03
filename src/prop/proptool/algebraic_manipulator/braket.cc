#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

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

  Total_Op_ = MT_map->at(multiop_name_);

  vector<char> projector_names; 
  for ( auto& tens : Total_Op_->sub_tensops() ) 
    if ( tens->is_projector() ) 
      projector_names.push_back(tens->name()[0]);

  print_vector( op_trans_list_ , "op_trans_list" ) ; cout << endl;
 
  auto trans_info = make_pair( op_trans_list_,  op_order_ );
  auto GGen = make_shared<GammaGeneratorRedux>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, G_to_A_map, factor_ );
  for ( auto range_map_it = Total_Op_->split_rxnges()->begin(); range_map_it !=Total_Op_->split_rxnges()->end(); range_map_it++ ){
//    if ( !range_map_it->second->is_sparse( op_state_ids_ ) ){    
     shared_ptr<Range_BlockX_Info>  trans_block = range_map_it->second->transform( make_shared<const vector<string>>(range_map_it->first), Total_Op_->idxs(), Total_Op_->aops(), op_order_, op_trans_list_ );
      
      GGen->add_gamma( trans_block , range_map_it->first );
      
      if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){
        if ( GGen->generic_reorderer( "normal order", false, false ) ) {
          if ( GGen->generic_reorderer( "alternating order", false, true ) ){  

//            cout << "We need these blocks : " ; cout.flush(); cout << " Total_Op_->sub_tensops().size() = " ; cout.flush(); cout << Total_Op_->sub_tensops().size() << endl; 
//            vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();
//            int qq = 0 ;
//            for ( auto& tens_block : *(range_map_it->second->range_blocks()) ){ 
//              cout << tens_block->orig_name()  << " " ; cout.flush(); cout << "sub_tensops[" << qq<< " ]->name() = "; cout.flush(); cout << sub_tensops[qq]->name() << endl;
//              MT_map->at( sub_tensops[qq++]->name()  )->add_required_block( tens_block->orig_name() );
//              required_blocks->emplace( tens_block->orig_name() );
//            }
//            cout << endl;
//          }
          cout << "got a thing" << endl;   
        }
      }
    }
  }

  print_gamma_Atensor_contractions( G_to_A_map, false );

  return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::print_gamma_Atensor_contractions(shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo> >>>> G_to_A_map,
		                                        bool has_orb_exc ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout <<  "BraKet<DataType>::print_gamma_Atensor_contractions()" << endl; 

  if ( !has_orb_exc ) { 
    cout << "no proj" << endl;
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

  } else {
    cout << "no proj" << endl;
    for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){

      cout << "====================================================" << endl;
      cout << map_it->first << endl;
      cout << "====================================================" << endl;

      if ( map_it->second->size() == 0 ) {

        cout << endl << "      - No Contributions! -" <<  endl << endl;

      } else {

      for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
        cout <<  A_map_it->first << "  "; cout.flush();
        string spacer = "  "; for ( int ss =0 ; ss != A_map_it->first.size() ; ss++ )  spacer += ' '; 
          for ( int qq = 0; qq !=A_map_it->second->aid_orders().size(); qq++) {
            cout << "{ " ; print_vector( *(A_map_it->second->aid_order(qq)), "" );
            cout << " : [ "; cout.flush();
            for ( int rr = 0; rr != A_map_it->second->aid_pid_orders(qq)->size(); rr++ ) {
              print_vector( *(A_map_it->second->pid_order(qq,rr)), "" ); cout.flush();
              cout << "*(" << A_map_it->second->factor(qq,rr).first << "," <<  A_map_it->second->factor(qq,rr).second << ") "; cout.flush();
            }
            cout << "] }" << endl;
            cout << spacer ; cout.flush();
          }
          cout << endl;
        }
      }
    }
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<double>;
//template class BraKet<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
