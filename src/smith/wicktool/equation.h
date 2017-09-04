#include <src/smith/wicktool/BraKet.h>
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/TensOp.h>
#include <src/smith/tensor.h>
 #ifndef __SRC_SMITH_EQUATION_H
 #define __SRC_SMITH_EQUATION_H

//#include "WickUtils.h"
//#include "BraKet.h"
//#include "TensOp.h"


template<class DType> 
class Equation {
      public:

      Equation(){};
      ~Equation(){};

      int nidxs = 8;
      int nact = 5;
      int norb = 5;
      bool spinfree = false;

      //Equation builder
      void equation_build(std::shared_ptr<std::vector< std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> >> BraKet_list);

      //Vector of BraKet terms which comprise the equation
      std::vector<std::shared_ptr<BraKet<DType>>> BraKet_Terms;
    
      // tensor data map
      std::shared_ptr< std::map< std::string, std::shared_ptr<DType> >>Tparts_map ;

      //contracted and uncontracted single tensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DType> > >>CTP_map    ;      

      // contracted and uncontracted multitensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DType> > >>CMTP_map   ;  

      // Vector containing the operations which need to be performed to calculate a given A-tensor
      std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>,std::string> >> ACompute_list ;

      // key : String identifying the Gamma or rdm      result : Vector of string vectors identifying A-Tensor contributions paired with corresponding ReIm factor 
      std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> CMTP_Eqn_Compute_List ;

      void Add_BraKet_Compute_Terms_CMTP(std::shared_ptr<BraKet<DType>> BK );

      void Initialize();

      void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> Tens_vec  );
      
      std::shared_ptr<TensOp<DType>> Build_TensOp(std::string op_name,
                                std::shared_ptr<DType> tensor_data, //needs to be templated, should be Bagel tensor
                                std::shared_ptr<std::vector<std::string>> op_idxs,
                                std::shared_ptr<std::vector<bool>> op_aops, 
                                std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                                std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> Symmetry_Funcs,
                                std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> Constraint_Funcs,
                                std::pair<double,double> factor, std::string Tsymmetry, bool hconj ) ;

};
#endif
