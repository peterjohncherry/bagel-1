#ifndef __SRC_SMITH_WICKTOOL_GAMMA_COMPUTER_H
#define __SRC_SMITH_WICKTOOL_GAMMA_COMPUTER_H

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/wicktool/equation_computer.h>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/util/f77.h>
#include <src/smith/wicktool/tensor_sorter.h>
#include <src/smith/wicktool/tensor_arithmetic.h>

namespace bagel {
namespace SMITH { 

namespace Gamma_Computer { 

class Gamma_Computer { 

    public: 
    Gamma_Computer( std::shared_ptr<Equation<double>> Eqn_info_in,
                    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in      );
    ~Gamma_Computer(){};
  
    int nelea_  ;
    int neleb_  ;
    int norb_   ;
    int nstate_ ;
    int maxtile ;
    int cimaxblock;
    std::shared_ptr<IndexRange> virt_   ;
    std::shared_ptr<IndexRange> active_ ;
    std::shared_ptr<IndexRange> closed_ ;

    std::shared_ptr<const Dvec> cc_; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo>>> GammaMap;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Gamma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Sigma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CIvec_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<const Determinants>>> Determinants_map;

    std::shared_ptr<Equation<double>> Eqn_info;

    std::shared_ptr<Tensor_Arithmetic::Tensor_Arithmetic<double>> Tensor_Calc;

    void Initialize_wfn_info( std::shared_ptr<Civec> civector, int state_num ) ;

    void get_gamma_tensor( int MM , int NN, std::string gamma_name);

    void compute_gammas_blocked(const int MM, const int NN, std::string gamma_name) ;

    void get_gamma_2idx(const int MM, const int NN, std::string gamma_name );

    std::unique_ptr<double[]>
    gamma_2idx_block(std::shared_ptr<const Civec> cbra, std::shared_ptr<const Civec> cket, std::pair<size_t,size_t> irange,
                     std::pair<size_t,size_t> jrange ) const;

    std::unique_ptr<double[]> 
    sigma_blocked(std::shared_ptr<const Civec> cvec, std::pair<size_t,size_t> irange, std::pair<size_t,size_t> jrange) const ;
   
    std::shared_ptr<Tensor_<double>> 
    convert_civec_to_tensor( std::shared_ptr<const Civec> civector, int state_num ) const ;

    void get_civector_indexranges(int nstates) ;

    void build_sigma_block( std::string sigma_name, std::vector<Index>& sigma_id_blocks, std::vector<int>& sigma_offsets, std::string Ket_name  ) const ;

    void build_sigma_2idx_tensor(std::string Bra_name, std::string Ket_name, std::shared_ptr<std::vector<std::string>> orb_ranges_str)  ;

    void build_gamma_2idx_tensor( int NN, int MM, int nelea, int neleb, int norb, std::string gamma_name ) ;

    /////////// Utility routines /////////////////////////
 
    std::shared_ptr<std::vector<IndexRange>> Get_Bagel_IndexRanges( std::shared_ptr<std::vector<std::string>> ranges_str);

    std::string get_sigma_name( std::string Bra_name, std::string Ket_name , std::shared_ptr<std::vector<std::string>>  orb_ranges, std::shared_ptr<std::vector<bool>> aops ) ;

    std::string get_civec_name(int state_num, int norb, int nalpha, int nbeta) const ;
 
    std::string get_det_name(std::shared_ptr<const Determinants> Detspace ) ;

    void tester();

};

}
}
}
#endif
