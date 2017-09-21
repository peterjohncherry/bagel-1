#ifndef __SRC_SMITH_BraKet_H
#define __SRC_SMITH_BraKet_H
 #include <src/smith/wicktool/WickUtils.h>
 #include <src/smith/wicktool/gamma.h>
 #include <src/smith/wicktool/gamma_generator.h>
 #include <src/smith/wicktool/TensOp.h>

// #include "WickUtils.h"
// #include "gamma.h"
// #include "gamma_generator.h"
// #include "TensOp.h"

template<class DType> 
class BraKet  {
      using pint_vec = std::vector<std::pair<int,int>>;
      using pstr_vec = std::vector<std::pair<std::string,std::string>>;
      
      using CombinedGammaMap = std::map<std::vector<std::string>, /*spins of gamma indexes*/ 
                                        std::tuple< std::shared_ptr<std::vector< std::shared_ptr< std::vector<std::pair<std::string, std::string>> >>>, /* delta indexes  */  
                                                    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string, std::string>> >>>, /* spins of delta indexes */  
                                                    std::shared_ptr<std::vector<int>>  /* signs from reordering */ > >;

      using RelCombinedGammaMap = std::map< std::pair<std::vector<std::string>, std::pair<int,int>>,  /*spins of gamma_ids and original spinsector*/ 
                                            std::tuple< std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* contractions from reordering to gamma */  
                                                        std::shared_ptr<std::vector<std::shared_ptr<pstr_vec>>>, /* spins of contractions from reordering to gamma */  
                                                        std::shared_ptr<std::vector<int>>  /* signs from reordering to gamma */ > >;


      public:
      BraKet();
      ~BraKet(){};

     int nact_orb;
     int nact_el;
     int nidxs;

     bool spinfree_;
     int spin_max ;
     int spin_min ;

     std::shared_ptr<CombinedGammaMap> GammaMap;

     std::shared_ptr<RelCombinedGammaMap> RelGammaMap;
  
     std::shared_ptr<std::vector<std::shared_ptr<TensOp<DType>>>> Sub_Ops;
     std::shared_ptr<MultiTensOp<DType>> Total_Op;

     std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors ;

     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths;
     
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> >> spin_paths_str;

     std::shared_ptr< std::map< std::vector<std::string>, std::shared_ptr<std::vector<std::pair< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> >>> >> BK_Compute_List_CMTP;

    //functions 
    void add_Op(std::string op_name,
                    std::shared_ptr<std::vector<std::string>> op_idxs,
                    std::shared_ptr<std::vector<bool>> op_aops, 
                    std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                    std::vector<std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> Symmetry_Funcs,
                    std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> Constraint_Funcs,
                    std::pair<double,double> factor, std::string Tsymmetry, bool hconj ) ;

     void Build_TotalOp();

     void initialize(int nel, int norb, int nidxs, bool spinfree );

     void set_spin_info();
     
     std::shared_ptr<std::vector<std::pair<int,int>>> spin_sector_gen(int nact_el, int nact_orb);

     std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> all_spin_combinations(int length);

     std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> sector_constrained_spin_combs(int nalpha, int nbeta);

     bool sector_constrained_spin_combs(int nalpha, int nbeta , std::shared_ptr<std::vector<std::string>> gamma_comb);

     bool sector_constrained_spin_combs(int nalpha, int nbeta , std::shared_ptr<std::vector<bool>> gamma_comb);

     void bfunc(int xx, std::shared_ptr<std::vector<bool>> xvec, int nel, int length, std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> all_sc);
   
     void Build_Gamma_WithSpin(std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<std::vector<std::string>> idxs);

     void Build_Gamma_SpinFree(std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<std::vector<std::string>> idxs);

     void Build_Gamma_SpinFree_New(std::shared_ptr<std::vector<bool>> aops, std::shared_ptr<std::vector<std::string>> idxs);

     void get_string_spin_paths(); //should not be needed...

     void Build_Tensor_Contraction_list_CMTP();

//   void Build_Tensor_Contraction_list_Rel();
};
#endif
