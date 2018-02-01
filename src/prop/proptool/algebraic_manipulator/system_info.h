 #ifndef __SRC_PROP_PROPTOOL_SYSTEM_INFO_H
 #define __SRC_PROP_PROPTOOL_SYSTEM_INFO_H
 #include <src/prop/proptool/algebraic_manipulator/expression.h>
 #include <src/prop/proptool/algebraic_manipulator/states_info.h>
 #include <src/prop/proptool/algebraic_manipulator/equation.h>
 #include <src/smith/tensor.h>
using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

template<class DataType> 
class System_Info {
      private :
        bool spinfree_ = false;
        std::shared_ptr<StatesInfo<DataType>> target_states_;

      public:

      std::vector<std::string> free;
      std::vector<std::string> not_core;
      std::vector<std::string> not_act;
      std::vector<std::string> not_virt;
      std::vector<std::string> core;
      std::vector<std::string> act;
      std::vector<std::string> virt;

      // key :    Name of BraKet
      // result : Vector of TensOps corresponding to BraKet
      std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DataType>>>>>> braket_map_;
         
      // key :    Name of Term
      // result : Vector of BraKets corresponding to Term
      std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>> term_braket_map_;

      // key :    Name of Expression
      // result : Vector of Terms corresponding to Expression
      std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_;

      // key : equation_name 
      // result : equation object 
      std::shared_ptr< std::map <std::string, std::shared_ptr<Expression<DataType>>>> expression_map;

      // key : expression name
      // result : expresion object 
      std::shared_ptr< std::map <std::string, std::shared_ptr<Equation_Base<DataType>>>> equation_map_;
    
      // key :    Name of uncontracted part of TensorOp.
      // result : Info for uncontracted part of TensorOp info.
      std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp::TensOp<DataType> > >> T_map    ;      

      // key :    Name of uncontracted part of MultiTensorOp.
      // result : Info for uncontracted part of MultiTensorOp info.
      std::shared_ptr< std::map< std::string, std::shared_ptr< MultiTensOp::MultiTensOp<DataType> > >> MT_map    ;      

      // key :    Name of contracted part of TensorOp (CTP)
      // result : Info for contracted part of TensorOp info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart<DataType> > >> CTP_map;      

      // key :    Name of contracted part of multitensorop (CMTP)
      // result : Info for contracted part of multitensorop info
      std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> > >> CMTP_map   ;  

      // key : name of ATensor (also name of CTP)     
      // result : contraction list list for calculating that ATensor
      std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base>>>>> ACompute_map;

      // ket : name of Gamma
      // result : GammaInfo
      std::shared_ptr< std::map <std::string, std::shared_ptr<GammaInfo>>> Gamma_map;

      /////CONSTRUCTOR and DESTRUCTOR/////
      System_Info(std::shared_ptr<StatesInfo<DataType>> target_states_, bool spinfree);
      ~System_Info(){};

      std::shared_ptr<TensOp::TensOp<DataType>> Initialize_Tensor_Op_Info( std::string op_name ) ;

      void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>> Tens_vec  );
      
      std::shared_ptr<TensOp::TensOp<DataType>> Build_TensOp( std::string op_name,
                                                              std::shared_ptr<std::vector<std::string>> op_idxs,
                                                              std::shared_ptr<std::vector<bool>> op_aops, 
                                                              std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                                                              std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> Symmetry_Funcs,
                                                              std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> Constraint_Funcs,
                                                              DataType factor, std::string Tsymmetry, bool hconj, int state_dependence ) ;
     
      void Set_BraKet_Ops(std::shared_ptr<std::vector<std::string>> Op_names, std::string term_name ) ;

      std::string Build_Expression( std::shared_ptr<std::vector<BraKet<DataType>>>  term_info_list  ); 
      std::string Build_Expression( std::string expression_name  ); 
    
      void create_equation( std::string name, std::string type, 
                            std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<BraKet<DataType>>>>>  term_braket_map,
                            std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map );

      std::string Get_BraKet_name( BraKet<DataType>& BraKet_info  ) ;

      int nalpha(int state_num) { return target_states_->nalpha( state_num ); };
      int nbeta(int state_num)  { return target_states_->nbeta( state_num );  };
      int nact(int state_num)   { return target_states_->nact( state_num );   };
      bool spinfree(){return spinfree_;}

      static std::string flip(std::string idx);
      static std::shared_ptr<std::vector<std::string>> ijkl_to_klij(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> ijkl_to_jilk(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> ijkl_to_lkji(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> ijkl_to_ijlk_block(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> ijkl_to_jikl_block(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> ijkl_to_jilk_block(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> bbbb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> bbaa_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> aabb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
      static std::shared_ptr<std::vector<std::string>> identity(std::shared_ptr<std::vector<std::string>> invec) ;
      static bool NotAllAct(std::shared_ptr<std::vector<std::string>> ranges);
      static bool always_true(std::shared_ptr<std::vector<std::string>> ranges);

      std::vector<std::tuple<std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> set_2el_symmfuncs();
      std::vector<std::tuple<std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> set_1el_symmfuncs();
      std::vector<std::tuple<std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> identity_only();

      struct compare_string_length {
        bool operator()(const std::string& first, const std::string& second) {
            return first.size() > second.size();
        }
      };

};

#endif
