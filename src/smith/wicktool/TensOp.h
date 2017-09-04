#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
#include <src/smith/wicktool/WickUtils.h>
#include <src/smith/wicktool/CtrTensOp.h>

//#include "WickUtils.h"
//#include "CtrTensOp.h"

using namespace WickUtils;

template<class DType>
class TensOp {
   private:
     std::string Tsymm_;
     std::string psymm_;
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;

   public:
    
     std::string name_;
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_;
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints_;  
     std::shared_ptr<std::vector<std::string>> idxs;
     std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges;
     std::shared_ptr<std::vector<bool>> aops;
     std::shared_ptr<DType> data;
     bool contracted;
     bool spinfree;

     std::shared_ptr<std::vector<int>> plus_ops;
     std::shared_ptr<std::vector<int>> kill_ops;
     std::shared_ptr< std::map< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > > orb_ranges;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DType>> > > CTP_map ;

     int factor;
     void get_ctrs_pos() ;
     
     //map from original index range, to ranges and factors used in calculation
     //tuple contained 1: bool is_this_range_unique?, 2: unique range which needs to be calculated, 3: indexes for contraction of this range,  4: factor from transformation 
     std::shared_ptr<std::map< std::shared_ptr<std::vector<std::string>>, 
                     std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges_;


     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    plus_combs;
     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    kill_combs;

   public: 
     TensOp() {};
     TensOp(std::string name,
            std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs, 
            std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints) : name_(name), symmfuncs_(symmfuncs), constraints_(constraints) {};

     ~TensOp(){};
   
     virtual void initialize(std::vector<std::string> orig_idxs, std::vector<std::vector<std::string>> orig_idx_ranges,
                             std::vector<bool> orig_aops, int orig_factor, std::string orig_Tsymm, std::string orig_psymm = "2el" );
     
     virtual
     std::shared_ptr<std::map< std::shared_ptr<std::vector<std::string>>, 
                     std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges(){ return all_ranges_; }

     void hconj();
  
     std::string name() {return name_ ;}   
  
     void generate_ranges( std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges, std::string symm_type);

     bool apply_symmetry( std::shared_ptr<std::vector<std::string>> ranges_1, std::shared_ptr<std::vector<std::string>> ranges_2  );

     bool satisfies_constraints( std::shared_ptr<std::vector<std::string>> ranges );

     void generate_all_contractions();
     
     virtual void get_ctrs_tens_ranges() ;
 
}; 

template<class DType>
class MultiTensOp : public TensOp<DType> {
   public: 
 
//     using TensOp<DType>::name_;
//     using TensOp<DType>::plus_ops;
//     using TensOp<DType>::kill_ops;
//     using TensOp<DType>::orb_ranges;
//     using TensOp<DType>::CTP_map ;
//     using TensOp<DType>::factor;
//     using TensOp<DType>::data; 
//     using TensOp<DType>::contracted;
//     using TensOp<DType>::spinfree;
//     using TensOp<DType>::aops;
//     using TensOp<DType>::idx_ranges;
//     using TensOp<DType>::idxs;




     std::string name_;
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_;
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints_;  
     std::shared_ptr<std::vector<std::string>> idxs;
     std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges;
     std::shared_ptr<std::vector<bool>> aops;
     std::shared_ptr<DType> data;
     bool contracted;
     bool spinfree;

     std::shared_ptr<std::vector<int>> plus_ops;
     std::shared_ptr<std::vector<int>> kill_ops;
     std::shared_ptr< std::map< std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > > orb_ranges;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DType>> > > CTP_map ;

     int factor;

     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    plus_combs;
     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    kill_combs;

   public: 
     MultiTensOp() {};
     MultiTensOp(std::string name, bool spinfree) : name_(name), spinfree_(spinfree) {}; // spin free construction
     MultiTensOp<DType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};

     ~MultiTensOp() {};





   private:
     //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
     std::shared_ptr<std::map< std::shared_ptr<std::vector<std::string>>,
                               std::tuple<bool, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> combined_ranges;

     //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
     std::shared_ptr<std::map< std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>>, 
                     std::tuple< std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> split_ranges;

   public: 

     std::shared_ptr<std::vector<int>> Tsizes ;
     std::shared_ptr<std::vector<int>> cmlsizevec ;
     std::shared_ptr<std::vector<std::string>> names;

     std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> split_idxs   ;
     std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> split_plus_ops;
     std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> split_kill_ops;

     std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors_;
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;
     bool spinfree_;

     //takes from pair of unc _range_ (defines gamma) and ctr_indexes (defines A contrib) 
     std::shared_ptr<std::map< std::tuple< std::vector<std::string>, pstr_vec, pstr_vec >, std::shared_ptr<std::vector<std::string>> >> CMTP_gamma_contribs;

     std::vector<std::shared_ptr<TensOp<DType>>> orig_tensors_;

     std::shared_ptr<std::map< std::string, std::shared_ptr<CtrMultiTensorPart<DType>> >> CMTP_map; 

    
    void initialize(std::vector<std::shared_ptr<TensOp<DType> >> orig_tensors);
    void generate_ranges();

    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors, std::shared_ptr<std::vector<std::string>> id_ranges );

    void get_ctrs_tens_ranges() override; 
   
    void print_gamma_contribs();
    
}; 
#endif
