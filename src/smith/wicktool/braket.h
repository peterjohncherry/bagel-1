#ifndef __SRC_SMITH_BraKet_H
#define __SRC_SMITH_BraKet_H
  #include <src/smith/wicktool/wickutils.h>
  #include <src/smith/wicktool/gamma_generator.h>
  #include <src/smith/wicktool/tensop.h>
  #include <src/smith/wicktool/spin_manager.h>
  #include <src/smith/wicktool/states_info.h>

 //#include "wickutils.h"
 //#include "gamma_generator.h"
 //#include "tensop.h"
 //#include "spin_manager.h"

template<class DType> 
class BraKet  {
      using pint_vec = std::vector<std::pair<int,int>>;
      using pstr_vec = std::vector<std::pair<std::string,std::string>>;

      public:
       BraKet(std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map_in,
              std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> GammaMap ,
              std::shared_ptr<StatesInfo<DType>> TargetStates  );
       ~BraKet(){};
      
       std::shared_ptr<StatesInfo<DType>> TargetStates;
       int nact_orb;
       int nact_el;
       int nidxs;
       
       bool spinfree;
       int spin_max ;
       int spin_min ;
       
       std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> Sub_Ops;
       
       std::shared_ptr<MultiTensOp<DType>> Total_Op;
       
       std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> BK_Compute_List_CMTP;
       
       //Takes each gamma name to a map containing the names of all A-tensors with which it must be contracted, and the relevant factors
       std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map;
       
       //Goes from a gamma name to the gamma info
       std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>> GammaMap;
       
       //functions 
       void add_Op(std::string op_name,
                   std::shared_ptr<std::vector<std::string>> op_idxs,
                   std::shared_ptr<std::vector<bool>> op_aops, 
                   std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                   std::vector<std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> Symmetry_Funcs,
                   std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> Constraint_Funcs,
                   std::pair<double,double> factor, std::string Tsymmetry, bool hconj ) ;
       
       void Build_TotalOp();
       
       void Build_Gamma_WithSpin(std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<std::vector<std::string>> idxs);
       
       void Build_Gamma_SpinFree(std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<std::vector<std::string>> idxs); 
       void Build_Tensor_Contraction_list_CMTP();

};
#endif
