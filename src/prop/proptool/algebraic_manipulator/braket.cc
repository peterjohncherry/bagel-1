#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_generator_redux.h>

using namespace std;
using namespace WickUtils;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
BraKet<DataType>::BraKet( std::vector<std::string>& op_list, std::vector<char>& op_trans_list,
                          DataType factor, int bra_num, int ket_num, 
                          std::shared_ptr<std::vector<std::vector<int>>> op_state_ids, std::string type) :
                          op_list_(op_list), op_trans_list_(op_trans_list), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
                          op_state_ids_(op_state_ids), type_(type), multiop_name_(std::accumulate(op_list_.begin(),
                          op_list_.end(), std::string(""))), proj_op_(false) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::BraKet" << endl;

  //get state name first
  string multiop_state_name = "";
  shared_ptr<vector<shared_ptr<vector<int>>>> state_id_list =  make_shared<vector<shared_ptr<vector<int>>>>();

  shared_ptr<vector<shared_ptr<Op_Info>>> multiop_info_list = make_shared<vector<shared_ptr<Op_Info>>>( op_list.size());
  vector<shared_ptr<Op_Info>>::iterator mil_it = multiop_info_list->begin();

  for ( int ii = 0 ; ii != op_list_.size(); ii++, mil_it++ ) {
    string op_state_name = "";
    op_state_name += op_list_[ii] ;
    if (op_state_ids_->at(ii).size() > 0 ) {
      op_state_name +=  "_{"; 
      for( int jj = 0; jj != op_state_ids_->at(ii).size(); jj++ ) 
        op_state_name += to_string(op_state_ids_->at(ii)[jj]); 
      op_state_name += "}"; 
    }

//    if (op_trans_list.size() > 0 ) {
//      op_state_name +=  "^{"; 
//      op_state_name += op_trans_list[ii]; 
//      op_state_name += "}"; 
//    }

    multiop_state_name += op_state_name;
    *mil_it = make_shared<Op_Info>( op_state_name, make_shared<vector<int>> (op_state_ids_->at(ii)), op_trans_list[ii] ); 
  }
  
  if (type_[0] == 'c' )// checking if derivative  
    name_ = "c_{I}"; 

  // Getting the BraKet name 
  name_ = "<" + std::to_string(bra_num)+ "| ";   name_ += multiop_state_name;   name_ += " |"+ std::to_string(ket_num) + ">";

  op_order_ = std::vector<int>(op_list.size());
  iota( op_order_.begin(), op_order_.end(), 0);
  sort( op_order_.begin(), op_order_.end(), [ &op_list ] ( int i1, int i2) { return (bool)( op_list[i1] < op_list[i2]); });  

  multiop_info_ = make_shared<MultiOp_Info>( multiop_state_name , multiop_info_list, op_order_ );  

  WickUtils::print_vector( op_list, "op_list" ); std::cout << std::endl; 
  WickUtils::print_vector( op_order_, "op_order" ); std::cout << std::endl; 
  
  projected_bra_ = false;
  projected_ket_ = false;
} 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following restructuing this class is starting to look more redundant, however I think it is still useful for
//merging, symmetry checking and sparsity. As well as controlling the reordering 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void BraKet<DataType>::generate_gamma_Atensor_contractions( shared_ptr<map<string,shared_ptr<TensOp_Base>>> MT_map,                
                                                            shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo> >>>> G_to_A_map,
                                                            shared_ptr<map<string, shared_ptr< GammaInfo >>> gamma_info_map,
                                                            shared_ptr<StatesInfo<DataType>> target_states,
                                                            shared_ptr<set<string>> required_blocks,
                                                            shared_ptr<map<string, shared_ptr<CtrTensorPart_Base>>> ctp_map  ) {  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::generate_gamma_Atensor_contractions : " << name_ << endl;

  Total_Op_ = MT_map->at( multiop_name_ );

  vector<char> projector_names;
  for ( auto& tens : Total_Op_->sub_tensops() )
    if ( tens->is_projector() )
      projector_names.push_back(tens->name()[0]);

//  auto trans_info = make_pair( op_trans_list_,  op_order_ );
//  vector<char>::const_iterator otl_it = op_trans_list_.begin();
//  for ( vector<int>::iterator oo_it = op_order_.begin();  oo_it != op_order_.end();  oo_it++, otl_it++ )
//    if ( *otl_it != '0' )
//      Total_Op_->sub_tensops()[*oo_it]->generate_transformed_ranges(*otl_it);

  shared_ptr<vector<bool>> trans_aops = Total_Op_->transform_aops( op_order_,  op_trans_list_ );

  print_vector(*trans_aops , " got transformed aops" ); cout <<endl;
 
  shared_ptr<GammaGeneratorRedux> GGen = make_shared<GammaGeneratorRedux>( target_states, bra_num_, ket_num_, Total_Op_, gamma_info_map, G_to_A_map, factor_ );

  for ( auto op_info : *(multiop_info_->op_info_vec_) ) { cout << "op_info->name_ = " <<  op_info->name_ << endl;} 

  Total_Op_->generate_ranges( multiop_info_ ); 

  auto split_ranges = Total_Op_->state_specific_split_ranges_->at( multiop_info_->name_ ); 

  for ( auto range_map_it = split_ranges->begin(); range_map_it != split_ranges->end(); range_map_it++ ){

    // get transformed range blocks;
    pair<double, double> block_factor = make_pair( factor_, factor_ );

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

  ctp_map->insert( Total_Op_->CTP_map()->begin(), Total_Op_->CTP_map()->end() );

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
