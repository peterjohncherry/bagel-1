#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BraKet_Base::BraKet_Base( std::shared_ptr<MultiOp_Info> multiop_info, 
                          std::pair<double, double> factor, int bra_num, int ket_num,  std::string type) :
                          multiop_info_(multiop_info), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
                          type_(type) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet_Base::BraKet_Base" << endl;

  if (type_[0] == 'c' ) {
    name_ = "c_{I}"; 
  } else {
    name_ = ""; 
  }
  name_ += "<" + std::to_string(bra_num)+ "| ";   name_ += multiop_info->name_;   name_ += " |"+ std::to_string(ket_num) + ">";
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
BraKet<DataType>::BraKet( std::shared_ptr<MultiOp_Info> multiop_info, 
                          std::pair<double, double> factor, int bra_num, int ket_num,  std::string type) :
                          BraKet_Base( multiop_info, factor, bra_num, ket_num, type) {} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<TensOp_Base>>> MT_map,
                                                            shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo<DataType>> >>>> G_to_A_map,
                                                            shared_ptr<map<string, shared_ptr< GammaInfo<DataType> >>> gamma_info_map,
                                                            shared_ptr<StatesInfo<DataType>> target_states,
                                                            shared_ptr<set<string>> required_blocks,
                                                            shared_ptr<map<string, shared_ptr<CtrTensorPart_Base>>> ctp_map  ){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::generate_gamma_Atensor_contractions : " << name_ << endl;

  Total_Op_ = MT_map->at( multiop_info_->op_name_ );

  vector<char> projector_names;
  for ( auto& tens : Total_Op_->sub_tensops() )
    if ( tens->is_projector() )
      projector_names.push_back(tens->name()[0]);

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( *(multiop_info_->op_order_),  *(multiop_info_->transformations_) );
 
  Total_Op_->generate_ranges( multiop_info_ );

  shared_ptr<GammaGeneratorRedux<DataType>> GGen = make_shared<GammaGeneratorRedux<DataType>>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, G_to_A_map, factor_ );

  auto split_ranges = Total_Op_->state_specific_split_ranges_->at( multiop_info_->name_ );

  for ( auto range_map_it = split_ranges->begin(); range_map_it != split_ranges->end(); range_map_it++ ){

    // TODO if is only here for non-relativistic case ; constraints should be specified in input
    if ( !(range_map_it->second->ci_sector_transition_ ) ) {

      GGen->add_gamma( range_map_it->second, trans_aops );

      if ( GGen->generic_reorderer( "anti-normal order", true, false ) ){ 
        if ( GGen->generic_reorderer( "normal order", false, false ) ) {
          if ( GGen->generic_reorderer( "alternating order", false, true ) ){
      
            cout << "We need these blocks : " ; cout.flush(); cout << " Total_Op_->sub_tensops().size() = " ; cout.flush(); cout << Total_Op_->sub_tensops().size() << endl;
            vector<shared_ptr<TensOp_Base>> sub_tensops = Total_Op_->sub_tensops();
            int qq = 0 ;
            for ( auto& tens_block : *(range_map_it->second->range_blocks()) ){
              cout << tens_block->name()  << " from" ; cout.flush(); cout << tens_block->full_op_name() <<  " is a required block " << endl;
              MT_map->at( sub_tensops[qq++]->name()  )->add_required_block( tens_block->name() );
              required_blocks->emplace( tens_block->name() );
            }
            cout << endl;
          }
        }
      }
    }
  }

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

  print_gamma_Atensor_contractions( G_to_A_map, false );

  return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::print_gamma_Atensor_contractions(shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo<DataType>> >>>> G_to_A_map,
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
template class BraKet<std::complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
