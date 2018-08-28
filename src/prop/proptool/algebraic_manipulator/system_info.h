#ifndef __SRC_PROP_PROPTOOL_SYSTEM_INFO_H
#define __SRC_PROP_PROPTOOL_SYSTEM_INFO_H
#include <src/prop/proptool/algebraic_manipulator/equation_linearRM.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

template<class DataType> 
class System_Info {
      private :
        bool spinfree_ = false;
        std::shared_ptr<StatesInfo<DataType>> states_info_;

      public:

        std::vector<std::string> free;
        std::vector<std::string> not_core;
        std::vector<std::string> not_act;
        std::vector<std::string> not_virt;
        std::vector<std::string> core;
        std::vector<std::string> act;
        std::vector<std::string> virt;
         
        // key :    Name of range
        // result : prime for this range
        std::shared_ptr< std::map < char, long unsigned int >>  range_prime_map_; 

        // key :    Name of BraKet
        // result : Vector of TensOps corresponding to BraKet
        std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DataType>>>>>> braket_map_;
           
        // key :    Name of Term
        // result : Vector of BraKets corresponding to Term
        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>> term_braket_map_;
        
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
        std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp::TensOp<DataType> > >> T_map_;      
        
        // key :    Name of uncontracted part of MultiTensorOp.
        // result : Info for uncontracted part of MultiTensorOp info.
        std::shared_ptr< std::map< std::string, std::shared_ptr<TensOp_Base> >> MT_map_;      
        
        // key :    Name of Contracted Tensor Op Block (CTP)  or  Contracted Multi Tens Op Block (CMTP) 
        // result : Info object for CTP or CMTP
        std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base > >> CTP_map_;      
        
        // key : name of ATensor (also name of CTP)     
        // result : contraction list list for calculating that ATensor
        std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base>>>>> ACompute_map;
        
        // ket : name of Gamma
        // result : GammaInfo_Base
        std::shared_ptr< std::map <std::string, std::shared_ptr<GammaInfo_Base>>> Gamma_map;
        
        /////CONSTRUCTOR and DESTRUCTOR/////
        System_Info(std::shared_ptr<StatesInfo<DataType>> states_info, bool spinfree);
        ~System_Info(){};
        
        void construct_equation_task_list( std::string equation_name );
        
        void Build_BraKet(std::shared_ptr<std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>> Tens_vec  );
        
        std::shared_ptr<TensOp::TensOp<DataType>> Build_TensOp( std::string op_name,
                                                                std::shared_ptr<std::vector<std::string>> op_idxs,
                                                                std::shared_ptr<std::vector<bool>> op_aops, 
                                                                std::shared_ptr<std::vector<std::vector<std::string>>> op_idx_ranges,
                                                                std::vector<std::shared_ptr<Transformation>>& symmfuncs,
                                                                std::vector<std::shared_ptr<Constraint>>& constraints,
                                                                std::pair<double,double> factor, std::string Tsymmetry, bool hconj, int state_dependence ) ;
       
        void  create_equation( std::string name, std::string type, 
                               std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>>  term_braket_map,
                               std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map, 
                               std::shared_ptr<std::map<std::pair< std::string, std::vector<std::pair<std::string, int>>>, 
                                                                   std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>> term_braket_map_state_spec, 
                               std::shared_ptr<std::map< std::pair<std::string, std::vector<std::pair<std::string, int>>>, 
                                                         std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map_state_spec ); 
        
        void create_equation( std::string name, std::string type, 
                              std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>>  term_braket_map,
                              std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType,std::string>>>>> expression_term_map );
        
        int nalpha(int state_num) { return states_info_->nalpha( state_num ); };
        int nbeta(int state_num)  { return states_info_->nbeta( state_num );  };
        int nact(int state_num)   { return states_info_->nact( state_num );   };
        bool spinfree() {return spinfree_;}
        
        std::shared_ptr< std::map <std::string, std::shared_ptr<std::vector<std::shared_ptr< TensOp::TensOp<DataType>>>>>> braket_map(){return braket_map_;}
           
        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::shared_ptr<BraKet_Base>>>>> term_braket_map(){return term_braket_map_;}
        
        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector<std::pair<DataType, std::string>>>>> expression_term_map(){ return expression_term_map_;}
        
        std::shared_ptr< std::map <std::string, std::shared_ptr<Equation_Base<DataType>>>> equation_map(){ return equation_map_;}
        
        std::shared_ptr< std::map< std::string, std::shared_ptr< TensOp_Base >>> MT_map() { return MT_map_; }     
        
        std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base > >> CTP_map() { return CTP_map_;   }   
        
};

#endif
