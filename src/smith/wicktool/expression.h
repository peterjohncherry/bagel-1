#ifndef __SRC_SMITH_EQUATION_H
#define __SRC_SMITH_EQUATION_H
#include <src/smith/wicktool/braket.h>
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/tensop.h>
#include <src/smith/wicktool/states_info.h>
#include <src/smith/tensor.h>

template<class DataType>
class Expression {
      public:

      Expression(){};
      Expression( std::shared_ptr<std::vector<std::pair<std::string, DataType>>> BraKet_list,
                  std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DataType>>>>>> BraKet_map,
                  std::shared_ptr<StatesInfo<DataType>> TargetsInfo );

      ~Expression(){};

      std::shared_ptr<StatesInfo<DataType>> target_states_;

      // Vector of BraKet terms which comprise the expression
      std::vector<std::shared_ptr<BraKet<DataType>>> BraKet_Terms;
    
      // Map containing original tensor operators info
      std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp::TensOp<DataType> > >> T_map;

      // Contracted and uncontracted single tensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DataType> > >> CTP_map;

      // Maps from BraKet name to list of operators for calculating that BraKet
      std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DataType>>>>>> BraKet_map;

      // contracted and uncontracted multitensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> > >> CMTP_map;

       // map from the name of a tensor, to the list of contractions which need to be performed to obtain it NEW
      std::shared_ptr<std::map<std::string,  std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >> ACompute_map;

      // key : String identifying the Gamma or rdm      result : Vector of string vectors identifying A-Tensor contributions paired with corresponding ReIm factor
      std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> CMTP_Eqn_Compute_List;

      //Takes each gamma name to a map containing the names of all A-tensors with which it must be contracted, and the relevant factors
      std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo > >>> G_to_A_map;
       
      std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo> > > GammaMap;

      void Initialize();

      void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>> Tens_vec, DataType factor  );
     
      void Get_CMTP_Compute_Terms();

};
#endif
