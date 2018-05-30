#ifndef __SRC_PROPTOOL
#define __SRC_PROPTOOL

#include <src/global.h>
#include <src/wfn/get_energy.h>
#include <src/wfn/reference.h>
#include <src/wfn/ciwfn.h>
#include <src/smith/indexrange.h>
#include <src/smith/tensor.h>
#include <src/smith/multitensor.h>
#include <src/prop/proptool/initialization/equation_init_linearRM.h>
#include <src/prop/proptool/initialization/op_bk_term_expr_init.h>
#include <src/prop/proptool/integrals/moint_computer.h>
#include <src/prop/proptool/task_translator/system_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/b_gamma_computer.h>


namespace bagel {
namespace PropTool { 

  class PropTool {

    // user input
    std::shared_ptr<const PTree> idata_;

    // information from reference wavefunction
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const CIWfn> ciwfn_;
    std::shared_ptr<const Dvec> civectors_;
    bool breit_;
    bool gaunt_;
    bool block_diag_fock_;
    bool spinfree_;
    int nstates_;
    std::vector<int> target_states_; //note this is total range of the states
    std::vector<int> all_states_;

    std::vector<std::vector<int>> degenerate_states_; // use as a hack for time reversal symmetry

    //range info 
    int nclosed_; 
    int ncore_; 
    int nact_; 
    int nvirt_; 
    int nocc_; 
    int nfrozenvirt_;
    size_t maxtile_;
    size_t cimaxtile_;

    //SMITH::Tensor range definitions
    std::shared_ptr<SMITH::IndexRange> closed_rng_; 
    std::shared_ptr<SMITH::IndexRange> active_rng_;  
    std::shared_ptr<SMITH::IndexRange> virtual_rng_;
    std::shared_ptr<SMITH::IndexRange> free_rng_  ; 
    std::shared_ptr<SMITH::IndexRange> not_closed_rng_ ; 
    std::shared_ptr<SMITH::IndexRange> not_active_rng_  ;
    std::shared_ptr<SMITH::IndexRange> not_virtual_rng_ ;

    std::shared_ptr<std::map< std::string, std::shared_ptr<SMITH::IndexRange>>> range_conversion_map_ ;
 
    //TODO come up with a better encoding
    std::vector<long unsigned int> range_primes_;
    std::vector<long unsigned int> ctr_primes_;
    std::shared_ptr<std::map< char, long unsigned int >> range_prime_map_;

  
    std::shared_ptr<System_Info<double>> sys_info_;
    std::shared_ptr<System_Computer::System_Computer<double>> system_computer_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression<double>>>> expression_map_;

    std::shared_ptr<std::map< std::string, std::shared_ptr<Term_Init> >> term_init_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Expression_Init> >> expression_init_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Equation_Init_Base> >> equation_init_map_;

    std::shared_ptr<StatesInfo<double>> targets_info_; //this contains all information about the states

    // initialiation info
    void set_primes();
    void read_input_and_initialize(); 
    void get_wavefunction_info();
    void set_target_state_info();
    void identify_degeneracies( const std::vector<double>& energies );
    void set_ao_range_info();
    void set_ci_range_info();

    std::shared_ptr<std::vector<SMITH::IndexRange>> convert_to_indexrange( std::shared_ptr<const std::vector<std::string>> range_block_str ) ;

    void get_terms_init( std::shared_ptr<const PTree> expression_inp ); 
    void get_equations_init( std::shared_ptr<const PTree> expression_init ); 
    void get_equation_init_LinearRM( std::shared_ptr<const PTree> equation_inp );
    void get_equation_init_Value( std::shared_ptr<const PTree> equation_inp );

    void get_new_ops_init( std::shared_ptr<const PTree> ops_def_tree ) ;
    void get_expression_variables( std::shared_ptr<const PTree> expr_def_tree ) ;

    void build_algebraic_task_lists( std::string  eqn_interdependence );
 
    std::shared_ptr<std::map< std::string , double >> inp_factor_map_;
    std::shared_ptr<std::map< std::string , std::shared_ptr<std::vector<double>> >> inp_indexed_factor_map_;
    std::shared_ptr<std::map< std::string , std::shared_ptr<std::vector<int>> >> inp_range_map_;

    std::vector<std::string> equation_execution_list_;

    void print_input_info();

    public: 

      PropTool(std::shared_ptr<const PTree> idata, std::shared_ptr<const Geometry> g, std::shared_ptr<const Reference> r);
      PropTool(); 
      ~PropTool(){};
    
     void compute() { std::cout << " not connected to anything yet" << std::endl;}; 
     void execute_compute_lists(); 
     int nclosed(){ return nclosed_;}     
     int nocc ()  { return nocc_;}     
     int nact ()  { return nact_;}     
     int ncore()  { return ncore_;}     
     int nvirt()  { return nvirt_;}     
     int nfrozenvirt(){ return nfrozenvirt_;}     
     bool block_diag_fock() { return block_diag_fock_; }     
 
};
};
};
#endif
