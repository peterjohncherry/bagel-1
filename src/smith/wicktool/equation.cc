#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation.h>

//#include "equation.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class DType>
void Equation<DType>::Initialize(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  CTP_map    = make_shared<map< string, shared_ptr<CtrTensorPart<DType>> >>();    
  CMTP_map   = make_shared<map< string, shared_ptr<CtrMultiTensorPart<DType>> >>(); 
  Tparts_map = make_shared<map< string, shared_ptr<DType> >>();                   
  
  CMTP_Eqn_Compute_List = make_shared<map< vector<string>, shared_ptr<vector<pair<shared_ptr<vector<string>>, pair<int,int> >>> >>();

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
shared_ptr<TensOp<DType>> Equation<DType>::Build_TensOp(string op_name,
                                          shared_ptr<DType> tensor_data, //needs to be templated, should be Bagel tensor
                                          shared_ptr<vector<string>> op_idxs,
                                          shared_ptr<vector<bool>> op_aops, 
                                          shared_ptr<vector<vector<string>>> op_idx_ranges,
                                          vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Symmetry_Funcs,
                                          vector<bool(*)(shared_ptr<vector<string>>)> Constraint_Funcs,
                                          pair<double,double> factor, string Tsymmetry, bool hconj ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp::BuildTensOp" <<   endl;

  if(hconj){
    reverse(op_idxs->begin(), op_idxs->end());  
    reverse(op_aops->begin(), op_aops->end());  
    reverse(op_idx_ranges->begin(), op_idx_ranges->end());  
    factor.second = -factor.second;// change this 
  }
  //NOTE: change to use proper factor
  int tmpfac =1;
  shared_ptr<TensOp<DType>>  New_Op = make_shared<TensOp<DType>>(op_name, Symmetry_Funcs, Constraint_Funcs);

  New_Op->data = tensor_data;
  New_Op->initialize(*op_idxs, *op_idx_ranges, *op_aops, tmpfac, Tsymmetry);
  New_Op->get_ctrs_tens_ranges();

  return New_Op;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adds terms associated with each gamma into the map
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void Equation<DType>::Add_BraKet_Compute_Terms_CMTP(shared_ptr<BraKet<DType>> BK  ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Add_BraKet_Compute_Terms_CMTP" << endl;  

  for (auto mapit = BK->GammaMap->begin(); mapit != BK->GammaMap->end(); mapit++) {
    for (int kk = 0 ; kk!= (get<0>(mapit->second))->size(); kk++){
      auto Acontrib_loc = BK->Total_Op->CMTP_gamma_contribs->find( tie(mapit->first, *(get<0>(mapit->second))->at(kk), *(get<1>(mapit->second))->at(kk) ));  
      if ( Acontrib_loc == BK->Total_Op->CMTP_gamma_contribs->end() ) 
        continue;

      for (auto CMTP_name : *Acontrib_loc->second){
        cout << endl<< CMTP_name <<endl; 
        CMTP_map->at(CMTP_name)->FullContract(CTP_map);
      }

      auto gamma_factor = make_pair( (get<2>(mapit->second))->at(kk),(get<2>(mapit->second))->at(kk) );
      auto BraKet_Contrib = make_pair(Acontrib_loc->second, gamma_factor);

      auto loc_in_map = CMTP_Eqn_Compute_List->find(mapit->first);  
      if ( loc_in_map == CMTP_Eqn_Compute_List->end() ) {
        auto init = make_shared<vector< pair<shared_ptr<vector<string>>,pair<int,int>> >>(1,BraKet_Contrib);
        CMTP_Eqn_Compute_List->emplace(mapit->first, init );
      } else {
        loc_in_map->second->push_back(BraKet_Contrib);
      }
    }
  }

  for (auto mapit = BK->GammaMap->begin(); mapit != BK->GammaMap->end(); mapit++) {
    print_vecX<string>(mapit->first , "gamma_spins");
    for (int kk = 0 ; kk!= (get<0>(mapit->second))->size(); kk++){
      print_pairvec<string>(*(get<0>(mapit->second))->at(kk), "gamma_ctrs");
      print_pairvec<string>(*(get<1>(mapit->second))->at(kk), "gamma_ctr_spins"); cout <<endl;
    }
  } 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void Equation<DType>::Build_BraKet(shared_ptr<vector<shared_ptr<TensOp<DType>>>> Tens_vec ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto New_Term = make_shared<BraKet<DType>>();
  New_Term->Sub_Ops = Tens_vec;
 
  New_Term->Build_TotalOp();
  New_Term->initialize(nact, norb, nidxs, spinfree );
//  New_Term->Total_Op->print_gamma_contribs();
  New_Term->Build_Gamma_SpinFree(New_Term->Total_Op->aops, New_Term->Total_Op->idxs); 

  CMTP_map->insert(New_Term->Total_Op->CMTP_map->begin(), New_Term->Total_Op->CMTP_map->end());
  BraKet_Terms.push_back(New_Term);   

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DType>
void  Equation<DType>::equation_build(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Into Equation" <<endl;

  Initialize();

  //////////////////Initialization section/////////////////////////
  ///////////This should be read in from an input file/////////////
  
  //spinfree orbital ranges
  vector<string> free     = {"cor", "act", "vir"};
  vector<string> not_core = {"act", "vir"};
  vector<string> not_act  = {"cor", "vir"};
  vector<string> not_virt = {"cor", "act"};
  vector<string> core = {"cor"};
  vector<string> act  = {"act"};
  vector<string> virt = {"vir"};
 
  ///////////////////////////////////X Tensor ////////////////////////////////////////
  string X_TimeSymm = "none";
  auto X_factor = make_pair(1.0,1.0);
  auto X_idxs = make_shared<vector<string>>(vector<string> {"X0", "X1", "X2", "X3"});
  auto X_aops = make_shared<vector<bool>>(vector<bool>  {true, true, false, false}); 
  auto X_idx_ranges = make_shared<vector<vector<string>>>( vector<vector<string>> { free, free, free, free }); 
  auto X_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)> X_constraints = { &always_true };
  shared_ptr<DType> X_data = make_shared<DType>();  

  auto XTens = Build_TensOp("X", X_data, X_idxs, X_aops, X_idx_ranges, X_symmfuncs, X_constraints, X_factor, X_TimeSymm, false ) ;
 
  /////////////////////////////////T Tensor ////////////////////////////////////////
  string T_TimeSymm = "none";
  auto T_factor = make_pair(1.0,1.0);
  auto T_idxs = make_shared<vector<string>>(vector<string>{"T0", "T1", "T2", "T3"}  );
  auto T_aops = make_shared<vector<bool>>  (vector<bool>  {true, true, false, false} ); 
  auto T_idx_ranges =  make_shared<vector<vector<string>>>( vector<vector<string>> { not_core, not_core, not_virt, not_virt });   
  auto T_symmfuncs = set_2el_symmfuncs();
  vector<bool(*)(shared_ptr<vector<string>>)> T_constraints = { &NotAllAct };
  shared_ptr<DType> T_data = make_shared<DType>();  

  auto TTens = Build_TensOp("T", T_data, T_idxs, T_aops, T_idx_ranges, T_symmfuncs, T_constraints, T_factor, T_TimeSymm, false ) ;
  ///////////////////////////////////////////////////////////////////////

  auto BraKet_Tensors = make_shared<vector<shared_ptr<TensOp<DType>>>>(vector<shared_ptr<TensOp<DType>>>{ XTens, TTens} );

  Build_BraKet( BraKet_Tensors );
  
  for (auto braket : BraKet_Terms){
    Add_BraKet_Compute_Terms_CMTP( braket );
  }

  cout << "Leaving Equation" <<endl;

  return ;
}

//int main() {
//
//   Equation<vector<double>> CIDeriv;
//
//   CIDeriv.equation_build();
//
//   return 1; 
//}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class Equation<std::vector<double>>;
template class Equation<bagel::SMITH::Tensor_<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
