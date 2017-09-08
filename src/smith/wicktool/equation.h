#include <src/smith/wicktool/BraKet.h>
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/TensOp.h>
#include <src/smith/tensor.h>
 #ifndef __SRC_SMITH_EQUATION_H
 #define __SRC_SMITH_EQUATION_H

//#include "WickUtils.h"
//#include "BraKet.h"
//#include "TensOp.h"

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

using CombinedGammaMap   =  std::map<std::vector<std::string>, /*spins of gamma indexes*/ 
                                     std::tuple< std::shared_ptr<std::vector< std::shared_ptr< std::vector<std::pair<std::string, std::string>> >>>, /* delta indexes  */  
                                                 std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string, std::string>> >>>, /* spins of delta indexes */  
                                                 std::shared_ptr<std::vector<int>>  /* signs from reordering */ > >;

using RelCombinedGammaMap = std::map< std::pair<std::vector<std::string>, std::pair<int,int>>,  /*spins of gamma_ids and original spinsector*/ 
                                      std::tuple< std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* contractions from reordering to gamma */  
                                                  std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* spins of contractions from reordering to gamma */  
                                                  std::shared_ptr<std::vector<int>>  /* signs from reordering to gamma */ > >;


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
    
      //Map containing original tensor operators info
      std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp<DType> > >> T_map    ;      

      //contracted and uncontracted single tensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DType> > >>CTP_map    ;      

      // contracted and uncontracted multitensor info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DType> > >>CMTP_map   ;  

       // map from the name of a tensor, to the list of contractions which need to be performed to obtain it
      std::shared_ptr<std::map<std::string,  std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>,std::string> >> >>  ACompute_map ;

      // key : String identifying the Gamma or rdm      result : Vector of string vectors identifying A-Tensor contributions paired with corresponding ReIm factor 
      std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> CMTP_Eqn_Compute_List ;

      void Add_BraKet_Compute_Terms_CMTP(std::shared_ptr<BraKet<DType>> BK );

      void Initialize();

      void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> Tens_vec  );
      
      std::shared_ptr<TensOp<DType>> Build_TensOp(std::string op_name,
                                std::shared_ptr<DType> tensor_data,
                                std::shared_ptr<std::vector<std::string>> op_idxs,
                                std::shared_ptr<std::vector<bool>> op_aops, 
                                std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                                std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> Symmetry_Funcs,
                                std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> Constraint_Funcs,
                                std::pair<double,double> factor, std::string Tsymmetry, bool hconj ) ;

};
#endif
