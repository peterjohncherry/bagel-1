#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/ctrtensop.h>

 //#include "wickutils.h"
 //#include "ctrtensop.h"

//using namespace WickUtils;

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

namespace TensOp_Prep {
template<class DataType>
class TensOp_Prep {
   public:
     std::string name_;
     std::string Tsymm_;
     std::vector<std::string> idxs_;
     std::vector<std::vector<std::string>> idx_ranges_;
     std::vector<bool> aops_;
     
     std::vector<int> plus_ops_;
     std::vector<int> kill_ops_;
     int num_idxs_;
     
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;
     
     DataType orig_factor_;
    
     std::shared_ptr< const std::vector<std::string>> idxs_const_ptr; 
     
     // unique_range_blocks tells you what parts of your tensor actually need to be calculated and stored.
     std::vector< std::shared_ptr< const std::vector<std::string>> > unique_range_blocks;
     
     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     std::map< const std::vector<std::string>, 
               std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > > all_ranges_;
     
     TensOp_Prep( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
                  std::vector<bool>& aops, DataType orig_factor,
                  std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
                  std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
                  std::string& Tsymm );
//     TensOp_Prep( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp_Prep<DataType>>> orig_tensors );
     TensOp_Prep(){};
     
     void generate_ranges();
     bool apply_symmetry( const std::vector<std::string>& ranges_1, const std::vector<std::string>& ranges_2  );
     bool satisfies_constraints( std::vector<std::string>& ranges );


};
}
namespace TensOp_General {

 template<class DataType>
 class TensOp_General {
    public:
      const std::string name_;
      const std::string Tsymm_;
      const std::vector<std::string> idxs_;
      const std::vector<bool> aops_;
      
      const std::vector<int> plus_ops_;
      const std::vector<int> kill_ops_;
      const int num_idxs_;
      const DataType factor_;
     
      // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
      const std::map< const std::vector<std::string>,
                      std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >> all_ranges_;
     
      // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
      const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks_;

      const std::shared_ptr<const std::vector<std::string>> idxs_ptr;
      const std::shared_ptr<const std::vector<bool>> aops_ptr;
      const std::shared_ptr<const std::vector<int>> plus_ops_ptr;
      const std::shared_ptr<const std::vector<int>> kill_ops_ptr;
      const std::shared_ptr<const std::map< const std::vector<std::string>,
                                       std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                   std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>  >> all_ranges_ptr;
      const std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks_ptr;
 
   public:
     TensOp_General(TensOp_Prep::TensOp_Prep<DataType>& TOp ) : name_(TOp.name_), Tsymm_(TOp.Tsymm_), idxs_(TOp.idxs_), aops_(TOp.aops_), plus_ops_(TOp.plus_ops_),
                                                   kill_ops_(TOp.kill_ops_), num_idxs_(TOp.num_idxs_), factor_(TOp.orig_factor_),all_ranges_(TOp.all_ranges_), 
                                                   unique_range_blocks_(TOp.unique_range_blocks)  ,
                                                   idxs_ptr(std::make_shared<const std::vector<std::string>>(idxs_)),
                                                   aops_ptr(std::make_shared<const std::vector<bool>>(aops_)),
                                                   plus_ops_ptr(std::make_shared<const std::vector<int>>(plus_ops_)),
                                                   kill_ops_ptr(std::make_shared<const std::vector<int>>(kill_ops_)), 
                                                   all_ranges_ptr(std::make_shared<const std::map< const std::vector<std::string>,
                                                                  std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                                  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>(all_ranges_)), 
                                                   unique_range_blocks_ptr(std::make_shared<const std::vector< std::shared_ptr< const std::vector<std::string>>>>(unique_range_blocks_)) 
                                                   {};
     ~TensOp_General(){};

   std::string const name(){ return name_;}
   std::shared_ptr<const std::vector<std::string>>  idxs(){ return idxs_ptr;}
   std::shared_ptr<const std::vector<bool>>  aops(){ return aops_ptr;}
   
   std::string const Tsymm(){ return Tsymm_; }

   int num_idxs() {return num_idxs_; }

   std::shared_ptr<const std::vector<int>> plus_ops(){ return plus_ops_ptr;}
   std::shared_ptr<const std::vector<int>> kill_ops(){ return kill_ops_ptr;}
   
   DataType factor() { return factor_; };
   std::shared_ptr<const std::map< const std::vector<std::string>,
                   std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                   all_ranges() {return  all_ranges_ptr; }
   // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
   const std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks(){ return  unique_range_blocks_ptr;}

};
}

namespace TensOp {
template<class DataType>
class TensOp {
   protected:
   /* const */std::string name_;
   /* const */std::string Tsymm_;
   /* const */std::string psymm_;
   /* const */std::shared_ptr<std::vector<std::string>> idxs_;
   /* const */std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges_;
   /* const */std::shared_ptr<std::vector<bool>> aops_;
               
   /* const */std::shared_ptr<std::vector<int>> plus_ops_;
   /* const */std::shared_ptr<std::vector<int>> kill_ops_;

   DataType factor_;
   int num_idxs_;

   private:
     std::shared_ptr<TensOp_General::TensOp_General<DataType>> tensop_orig_;
     std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;

     //map from original index range, to ranges and factors used in calculation
     //tuple contained 1: bool is_this_range_unique?, 2: unique range which needs to be calculated, 3: indexes for contraction of this range,  4: factor from transformation     
   public:
     TensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp<DataType> >> orig_tensors); 
     TensOp(std::shared_ptr<TensOp_General::TensOp_General<DataType> > tensop_orig  ): tensop_orig_(tensop_orig)  {} ; 
     ~TensOp(){};
 
     const std::shared_ptr<const std::map< const std::vector<std::string>,
                      std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                      all_ranges() const  {return tensop_orig_->all_ranges();} ;

     int num_idxs() {return tensop_orig_->num_idxs_; }
     int factor() { return tensop_orig_->factor_; };

     std::string const name(){ return tensop_orig_->name();}
     std::string const Tsymm(){ return tensop_orig_->Tsymm_; }

     std::shared_ptr<const std::vector<bool>>  aops(){ return tensop_orig_->aops_ptr;}
     std::shared_ptr<const std::vector<std::string>>  idxs(){ return tensop_orig_->idxs_ptr;}
     std::shared_ptr<const std::vector<int>> plus_ops(){ return tensop_orig_->plus_ops_ptr;}
     std::shared_ptr<const std::vector<int>> kill_ops(){ return tensop_orig_->kill_ops_ptr;}
     
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     
     virtual void get_ctrs_tens_ranges() ;
     
}; 
}
namespace MultiTensOp {
template<class DataType>
class MultiTensOp : public TensOp::TensOp<DataType> {
  private: 
    // should be a better way; could define both MultiTensOp and TensOp as derived from non-templated class, as none of these depend on DataType
    using TensOp::TensOp<DataType>::name_;
    using TensOp::TensOp<DataType>::Tsymm_;
    using TensOp::TensOp<DataType>::psymm_;
    using TensOp::TensOp<DataType>::idxs_;
    using TensOp::TensOp<DataType>::idx_ranges_;
    using TensOp::TensOp<DataType>::aops_;
    
    using TensOp::TensOp<DataType>::plus_ops_;
    using TensOp::TensOp<DataType>::kill_ops_;
    using TensOp::TensOp<DataType>::num_idxs_;

    std::vector<std::shared_ptr<TensOp::TensOp<DataType>>> orig_tensors_;
    std::shared_ptr<std::vector<int>> Tsizes_ ;
    std::shared_ptr<std::vector<int>> cmlsizevec_ ;
    std::shared_ptr<std::vector<std::string>> names_;
 
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> split_idxs_;

    //cout maps from state and spin-sector to list of possible Aops
    std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > ,
                              std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;

    std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors_;
    bool spinfree_;
  
    int num_tensors_;

    std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;

  private:
    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
    std::shared_ptr< std::map< std::vector<std::vector<std::string>>, 
                               std::tuple<std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>>>>> split_ranges;

  public: 

    //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
    std::shared_ptr<std::map< std::vector<std::string>,
                              std::tuple<bool, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> combined_ranges;

    //takes from pair of unc _range_ (defines gamma) and ctr_indexes (defines A contrib) 
    std::shared_ptr<std::map< std::tuple< std::vector<std::string>, pstr_vec, pstr_vec >, std::shared_ptr<std::vector<std::string>> >> CMTP_gamma_contribs;
    std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
    std::shared_ptr< std::map< std::string, std::shared_ptr<CtrMultiTensorPart<DataType>> >> CMTP_map; 

    MultiTensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType> >> orig_tensors);
    ~MultiTensOp() {};
   
    std::string name(){ return name_;}
    std::shared_ptr<std::vector<std::string>> idxs(){ return idxs_;}
    std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
    std::shared_ptr<std::vector<bool>> aops(){ return aops_;}
    
    std::string Tsymm() {return Tsymm_; }
    std::string psymm(){ return psymm_; }
    
    std::shared_ptr<std::vector<int>> plus_ops(){ return plus_ops_;}
    std::shared_ptr<std::vector<int>> kill_ops(){ return kill_ops_;}
    
    std::vector<std::shared_ptr<TensOp::TensOp<DataType>>> orig_tensors(){ return orig_tensors_ ;}
    std::shared_ptr<std::vector<int>> Tsizes() { return Tsizes_ ;}
    std::shared_ptr<std::vector<int>> cmlsizevec() { return cmlsizevec_ ; }
    std::shared_ptr<std::vector<std::string>> names() { return names_; }
    
    void generate_ranges();
    
    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors, std::shared_ptr<std::vector<std::string>> id_ranges );
    
    void get_ctrs_tens_ranges() override; 
    
    void print_gamma_contribs();

}; 
}
#endif
//     With spin construction 
//
//     MultiTensOp::MultiTensOp<DataType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
//                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
//                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};
