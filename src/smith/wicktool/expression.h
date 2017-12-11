 #ifndef __SRC_SMITH_EQUATION_H
 #define __SRC_SMITH_EQUATION_H
 #include <src/smith/wicktool/braket.h>
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/tensop.h>
 #include <src/smith/wicktool/states_info.h>
 #include <src/smith/tensor.h>

 //#include "wickutils.h"
 //#include "braket.h"
 //#include "tensop.h"
 //
using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

template<class DType> 
class Expression {
      public:

      Expression(){};
      Expression( std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::shared_ptr< TensOp<DType>>>>>> BraKet_list,
                std::shared_ptr<StatesInfo<DType>> TargetsInfo);

      ~Expression(){};

      std::shared_ptr<StatesInfo<DType>> TargetStates;

      //Vector of BraKet terms which comprise the expression
      std::vector<std::shared_ptr<BraKet<DType>>> BraKet_Terms;
    
      //Map containing original tensor operators info
      std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp<DType> > >> T_map    ;      

      //contracted and uncontracted single tensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DType> > >>CTP_map    ;      

      // contracted and uncontracted multitensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DType> > >>CMTP_map   ;  

       // map from the name of a tensor, to the list of contractions which need to be performed to obtain it NEW
      std::shared_ptr<std::map<std::string,  std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> >>  ACompute_map ;

      // key : String identifying the Gamma or rdm      result : Vector of string vectors identifying A-Tensor contributions paired with corresponding ReIm factor 
      std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> CMTP_Eqn_Compute_List ;

      //Takes each gamma name to a map containing the names of all A-tensors with which it must be contracted, and the relevant factors
      std::shared_ptr<std::map<std::string, std::shared_ptr< std::map<std::string, AContribInfo > >>> G_to_A_map;
       
      std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo> > > GammaMap; 

      void Initialize();

      void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> Tens_vec  );
     
      void Get_CMTP_Compute_Terms();

};
#endif
