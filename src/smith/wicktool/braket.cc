#include <bagel_config.h>
#ifdef COMPILE_SMITH
  #include <src/smith/wicktool/braket.h>
 //#include "braket.h"

using namespace std;
using namespace WickUtils;


using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;
      
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
BraKet<DataType>::BraKet(shared_ptr<map<string,shared_ptr<map<string,AContribInfo>>>> G_to_A_map_in,
                      shared_ptr<map<string, shared_ptr< GammaInfo >>> GammaMap_in,
                      shared_ptr<StatesInfo<DataType>> TargetStates_in  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Sub_Ops = make_shared<vector<shared_ptr<TensOp::TensOp<DataType>>>>(0);
  G_to_A_map = G_to_A_map_in;
  GammaMap = GammaMap_in;
  TargetStates = TargetStates_in;

} 

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
 
  int Ket_num = 0; 
  int Bra_num = 0; 
 
  shared_ptr<vector<string>> idxs_buff  = make_shared<vector<string>>(*idxs );
  shared_ptr<vector<bool>> aops_buff  = make_shared<vector<bool>>(*aops );
  print_vector( *idxs_buff , "idxs_buff" ) ; 
  print_vector( *aops_buff , "    aops_buff" ) ; cout << endl;
  cout << "Total_Op->all_ranges()->size() = " << Total_Op->all_ranges()->size() << endl;
  for (auto range_map_it = Total_Op->all_ranges()->begin() ;  range_map_it !=Total_Op->all_ranges()->end(); range_map_it++){
    if ( range_map_it->second->survives() ) { 
      shared_ptr<GammaGenerator>  GGen = make_shared<GammaGenerator>(TargetStates, Bra_num, Ket_num, aops_buff, idxs_buff, GammaMap, G_to_A_map); 
      GGen->add_gamma(make_shared<vector<string>>(range_map_it->first), 1) ;
      GGen->norm_order();
      GGen->optimized_alt_order();
    }
  }
 
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){

    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;

    if ( map_it->second->size() == 0 ) { 
       
      cout << "No Contributions!" <<  endl;   
          
    } else {

      for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
        cout <<  A_map_it->first << "  "; cout.flush();
        for ( int qq = 0; qq !=A_map_it->second.id_orders.size(); qq++) {
          cout << "[ " ; 
          print_vector( A_map_it->second.id_order(qq), " " ); 
          cout << " (" << A_map_it->second.factor(qq).first <<  "," <<  A_map_it->second.factor(qq).first<<  ") ]" ; cout <<endl;
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
