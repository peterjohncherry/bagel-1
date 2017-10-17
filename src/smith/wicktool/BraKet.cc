#include <bagel_config.h>
#ifdef COMPILE_SMITH
  #include <src/smith/wicktool/BraKet.h>
 //#include "BraKet.h"

using namespace std;
using namespace WickUtils;

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;
      
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
BraKet<DType>::BraKet(std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< std::unordered_map<std::string, std::pair<int,int> > >>> G_to_A_map_in,
                      std::shared_ptr<std::unordered_map<std::string, std::shared_ptr< GammaInfo >>> GammaMap_in ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Sub_Ops = make_shared<vector<shared_ptr<TensOp<DType>>>>(0);
  G_to_A_map = G_to_A_map_in;
  GammaMap = GammaMap_in;
} 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::add_Op(string op_name,
                    shared_ptr<vector<string>> op_idxs,
                    shared_ptr<vector<bool>> op_aops, 
                    shared_ptr<vector<vector<string>>> op_idx_ranges,
                    vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                    vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                    pair<double,double> factor, string Tsymmetry, bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  if(hconj){
    reverse(op_idxs->begin(), op_idxs->end());  
    reverse(op_aops->begin(), op_aops->end());  
    reverse(op_idx_ranges->begin(), op_idx_ranges->end());  
    factor.second = -factor.second;// change this 
  }

  //NOTE: change to use proper factor
  int tmpfac =1;
  shared_ptr<TensOp<DType>>  New_Op = make_shared<TensOp<DType>>(op_name, Symmetry_Funcs, Constraint_Funcs);

  New_Op->initialize(*op_idxs, *op_idx_ranges, *op_aops, tmpfac, Tsymmetry);
  New_Op->get_ctrs_tens_ranges();

  Sub_Ops->push_back(New_Op);
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_TotalOp(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "BraKet::Build_TotalOp" << endl;
  string MT_name = "";
  for (auto Tens : *Sub_Ops)
    MT_name += Tens->name_;

  Total_Op = make_shared<MultiTensOp<DType>>( MT_name , true);
  Total_Op->initialize( *Sub_Ops );
  Total_Op->get_ctrs_tens_ranges() ;
  cout << "Built_total_op" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Gamma_SpinFree(shared_ptr<vector<bool>> aops, shared_ptr<vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_SpinFree_New" << endl;
  cout << "aops = " ; for (bool aop : *aops) { cout << aop << " " ; } cout << endl;
 
  auto aops_buff  = make_shared<vector<bool>>(*aops );
  for (auto range_map_it = Total_Op->combined_ranges->begin() ;  range_map_it !=Total_Op->combined_ranges->end(); range_map_it++){
    auto GGen = make_shared<GammaGenerator>(aops_buff, idxs, GammaMap, G_to_A_map); 
    GGen->add_gamma(range_map_it->first, 1) ;
    GGen->norm_order();
    GGen->alt_order();
  }
 
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){
    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;
    for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
      cout <<  A_map_it->first  << "  (" << A_map_it->second.first  << "," <<  A_map_it->second.second << ")" << endl;
    }
  }

  return; 
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Gamma_WithSpin(shared_ptr<vector<bool>> aops, shared_ptr<vector<string>> idxs){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Build_Gamma_WithSpin" << endl;
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pairs up each gamma term (defined by a spin pathway,and spin sector)  with relevant CtrMultiTensorParts .
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void BraKet<DType>::Build_Tensor_Contraction_list_CMTP(){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class BraKet<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
