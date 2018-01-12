#ifndef __SRC_SMITH_BraKet_H
#define __SRC_SMITH_BraKet_H
  #include <src/smith/wicktool/wickutils.h>
  #include <src/smith/wicktool/gamma_generator.h>
  #include <src/smith/wicktool/tensop.h>
  #include <src/smith/wicktool/spin_manager.h>
  #include <src/smith/wicktool/states_info.h>

template<typename DataType> 
class BraKet{

      public:
       
        //Takes each gamma name to a map containing the names of all A-tensors with which it must be contracted, and the relevant factors
        std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map;
       
        //Goes from a gamma name to the gamma info
        std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo>>> GammaMap;
       
        std::shared_ptr<StatesInfo<DataType>> TargetStates;

        //factor; this is specific to the expression to which this BraKet object belongs
        DataType factor_;

        BraKet(std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo >>>> G_to_A_map_in,
               std::shared_ptr<std::map<std::string, std::shared_ptr< GammaInfo >>> GammaMap_in ,
               std::shared_ptr<StatesInfo<DataType>> TargetStates_in, DataType factor ) : 
               G_to_A_map(G_to_A_map_in), GammaMap(GammaMap_in), TargetStates(TargetStates_in), factor_(factor)
               {
                 Sub_Ops = std::make_shared<std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>>(0);
               }
        
        ~BraKet(){};
       
        int nact_orb;
        int nact_el;
        int nidxs;
        
        bool spinfree;
        int spin_max ;
        int spin_min ;
        
        std::shared_ptr<std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>> Sub_Ops;
        
        std::shared_ptr<MultiTensOp::MultiTensOp<DataType>> Total_Op;
        
        std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> BK_Compute_List_CMTP;
       
        void Build_TotalOp();
        
        void Build_Gamma_WithSpin(std::shared_ptr<const std::vector<bool>> aops, std::shared_ptr<const std::vector<std::string>> idxs);

        void Build_Gamma_SpinFree(std::shared_ptr<const std::vector<bool>> aops, std::shared_ptr<const std::vector<std::string>> idxs); 

        void Build_Tensor_Contraction_list_CMTP();

};
#endif
