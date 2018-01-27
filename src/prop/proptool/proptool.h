#ifndef __SRC_PROPTOOL
#define __SRC_PROPTOOL

#include <src/global.h>
#include <src/wfn/get_energy.h>
#include <src/wfn/reference.h>
#include <src/wfn/ciwfn.h>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/multitensor.h>
#include <src/smith/wicktool/system_info.h>
#include <src/smith/wicktool/expression_computer.h>
#include <src/smith/wicktool/expression.h>
#include <src/smith/wicktool/term.h>
#include <src/prop/proptool/moint_init.h>
#include <src/prop/proptool/moint.h>
#include <src/prop/proptool/moint_computer.h>
#include <src/prop/proptool/equation.h>


namespace bagel {
namespace PropTool { 

  class PropTool {

    std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const CIWfn> ciwfn_;
    std::shared_ptr<const Dvec> civectors_;
 
    int nclosed_; 
    int ncore_; 
    int nact_; 
    int nvirt_; 
    int nocc_; 
    int nfrozenvirt_;
  
    bool breit_;
    bool gaunt_;
    bool block_diag_fock_;
   
    size_t maxtile_;
    size_t cimaxtile_;
 
    std::string method_;

    std::shared_ptr<System_Info<double>> sys_info_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression<double>>>> expression_map_;
    std::shared_ptr<SMITH::Expression_Computer::Expression_Computer<double>> expression_machine_;
    std::shared_ptr<std::map< std::string, double>> scalar_results_map_;

    std::shared_ptr<std::map< std::string, std::shared_ptr<Term_Init> >> term_init_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression_Init> >> expression_init_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Equation_Init_Base> >> equation_init_map_;

    std::vector<int> target_states_;
    std::vector<int> all_states_;
    std::shared_ptr<StatesInfo<double>> targets_info_; 

    void set_target_state_info();
    void set_ao_range_info();
    void set_ci_range_info();

    void build_op_tensors( std::vector<std::string>& expression_list ) ;
    std::shared_ptr<std::vector<SMITH::IndexRange>> convert_to_indexrange( std::shared_ptr<const std::vector<std::string>> range_block_str ) ;

    void get_terms_init( std::shared_ptr<const PTree> expression_inp ); 
    void get_equations_init( std::shared_ptr<const PTree> expression_init ); 
    void get_equation_init_LinearRM( std::shared_ptr<const PTree> equation_inp );

    void get_new_ops_init( std::shared_ptr<const PTree> ops_def_tree ) ;
    void get_expression_variables( std::shared_ptr<const PTree> expr_def_tree ) ;

    std::shared_ptr<std::map< std::string , double >> inp_factor_map_;
    std::shared_ptr<std::map< std::string , std::shared_ptr<std::vector<double>> >> inp_indexed_factor_map_;
    std::shared_ptr<std::map< std::string , std::shared_ptr<std::vector<int>> >> inp_range_map_;

    public: 

      PropTool(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r);
      PropTool(); 
      ~PropTool(){};
    
      void compute() { std::cout << " not connected to anything yet" << std::endl;}; 

      std::shared_ptr<SMITH::IndexRange> closed_rng_; 
      std::shared_ptr<SMITH::IndexRange> active_rng_;  
      std::shared_ptr<SMITH::IndexRange> virtual_rng_;
      std::shared_ptr<SMITH::IndexRange> free_rng_  ; 
      std::shared_ptr<SMITH::IndexRange> not_closed_rng_ ; 
      std::shared_ptr<SMITH::IndexRange> not_active_rng_  ;
      std::shared_ptr<SMITH::IndexRange> not_virtual_rng_ ;
      std::vector<SMITH::IndexRange> pt2_ranges_;
      std::vector<SMITH::IndexRange> pt2_ranges_herm_conj_;
      
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> sigma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> civec_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> gamma_data_map_;
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::Tensor_<double>>>> tensop_data_map_;
      
      std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_ ;
       
      std::vector<std::shared_ptr<SMITH::MultiTensor_<double>>> T2_all;
      std::vector<std::shared_ptr<SMITH::MultiTensor_<double>>> lambda_all;
      std::shared_ptr<SMITH::Tensor_<double>> F_1el_all;
      std::shared_ptr<SMITH::Tensor_<double>> H_1el_all;
      std::shared_ptr<SMITH::Tensor_<double>> H_2el_all;// only {occ, virt, occ, virt});
      std::shared_ptr<SMITH::Tensor_<double>> v2_; 

     //std::shared_ptr<std::vector<std::string>> identity( std::shared_ptr<std::vector<std::string>> invec ) { return invec; }

     int nclosed(){ return nclosed_;}     
     int nocc ()  { return nocc_;}     
     int nact ()  { return nact_;}     
     int ncore()  { return ncore_;}     
     int nvirt()  { return nvirt_;}     
     int nfrozenvirt(){ return nfrozenvirt_;}     
     bool block_diag_fock() { return block_diag_fock_; }     
 
     std::string method(){ return method_; }
     std::shared_ptr<const Geometry> geom(){ return geom_;} ;
     std::shared_ptr<const Hcore> hcore(){ return ref_->hcore();}
     std::shared_ptr<const RDM<1>> rdm1_av(){ return  ref_->rdm1_av();}
 
};
};
};
#endif
