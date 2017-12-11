#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/ctrtensop.h>

 //#include "wickutils.h"
 //#include "ctrtensop.h"

using namespace WickUtils;

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
     std::shared_ptr< std::map< const std::vector<std::string>, 
                      std::tuple<const bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, const std::pair<int,int> > >> all_ranges_;
     
     TensOp_Prep( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
                  std::vector<bool>& aops, DataType orig_factor,
                  std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
                  std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
                  std::string& Tsymm );
     TensOp_Prep(){};
     
     void generate_ranges();
     bool apply_symmetry( const std::vector<std::string>& ranges_1, const std::vector<std::string>& ranges_2  );
     bool satisfies_constraints( std::vector<std::string>& ranges );


};

 template<class DataType>
 class TensOp_General {
    private:
      const std::string name_;
      const std::string Tsymm_;
      const std::vector<std::string> idxs_;
      const std::vector<bool> aops_;
      
      const std::vector<int> plus_ops_;
      const std::vector<int> kill_ops_;
      const int num_idxs_;
    
     // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
     const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_ranges_;
     
     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     const std::map< const std::vector<std::string>,
                     std::tuple<const bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, const std::pair<int,int>  >> all_ranges_;
    
     TensOp_General(TensOp_Prep<DataType> TOp ) : name_(TOp.name_), Tsymm_(TOp.Tsymm_), idxs_(TOp.idxs_), aops_(TOp.aops_), plus_ops_(*TOp.plus_ops_), kill_ops_(*TOp.kill_ops_), 
                                                  num_idxs_(TOp.num_idxs_), all_ranges_(TOp.all_ranges_) {};
     ~TensOp_General(){};
 };



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
   int num_idxs_;

   private:
     std::shared_ptr<const TensOp_General<DataType>> tensop_orig_;
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_;
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints_;  
  
   /* const*/ DataType factor_;

     std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;
     std::shared_ptr< std::map< std::vector<std::string>, 
                      std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges_;
     
     void get_ctrs_pos() ;

     //map from original index range, to ranges and factors used in calculation
     //tuple contained 1: bool is_this_range_unique?, 2: unique range which needs to be calculated, 3: indexes for contraction of this range,  4: factor from transformation     
      

   public:
     TensOp( std::string name,
             std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs, 
             std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints,
             std::vector<std::string> orig_idxs, std::vector<std::vector<std::string>> orig_idx_ranges,                      
             std::vector<bool> orig_aops, int orig_factor, std::string orig_Tsymm, std::string orig_psymm = "2el");
     TensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp<DataType> >> orig_tensors); 
     TensOp(std::shared_ptr< const TensOp_General<DataType> > tensop_orig  ) : tensop_orig_(tensop_orig) {}; 
     ~TensOp(){};
  
     std::string name(){ return name_;}
     std::shared_ptr<std::vector<std::string>> idxs(){ return idxs_;}
     std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
     std::shared_ptr<std::vector<bool>> aops(){ return aops_;}
     
     std::string Tsymm(){ return Tsymm_; }
     std::string psymm(){ return psymm_; }

     int num_idxs(){return num_idxs_; }

     std::shared_ptr<std::vector<int>> plus_ops(){ return plus_ops_;}
     std::shared_ptr<std::vector<int>> kill_ops(){ return kill_ops_;}
     
     int factor() { return factor_; };
     
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> >,
                               std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >>spin_paths(){ return spin_paths_;}
     
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs()
     { return symmfuncs_;}
     
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints(){ return constraints_;}
     
     bool contracted;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     
     virtual // change the mapped type so is purpose built class; want to avoid tuples for the sake of clarity.
     std::shared_ptr<std::map< std::vector<std::string>, 
                     std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges(){ return all_ranges_; }
     
     //should be set up to generate hermitian conjugated operator TODO check to see if this is working
     void hconj();
     
     void generate_ranges();
     virtual void get_ctrs_tens_ranges() ;
     
     bool apply_symmetry( std::shared_ptr<std::vector<std::string>> ranges_1, std::shared_ptr<std::vector<std::string>> ranges_2  );
     bool apply_symmetry( std::vector<std::string> ranges_1, std::vector<std::string> ranges_2  );
     
     bool satisfies_constraints( std::shared_ptr<std::vector<std::string>> ranges );
     bool satisfies_constraints( std::vector<std::string> ranges );
     
     void generate_all_contractions();
 
}; 

template<class DataType>
class MultiTensOp : public TensOp<DataType> {
  private: 
    // should be a better way; could define both MultiTensOp and TensOp as derived from non-templated class, as none of these depend on DataType
    using TensOp<DataType>::name_;
    using TensOp<DataType>::Tsymm_;
    using TensOp<DataType>::psymm_;
    using TensOp<DataType>::idxs_;
    using TensOp<DataType>::idx_ranges_;
    using TensOp<DataType>::aops_;
    
    using TensOp<DataType>::plus_ops_;
    using TensOp<DataType>::kill_ops_;
    using TensOp<DataType>::num_idxs_;

    std::vector<std::shared_ptr<TensOp<DataType>>> orig_tensors_;
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

    MultiTensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp<DataType> >> orig_tensors);
    ~MultiTensOp() {};
   
    std::string name(){ return name_;}
    std::shared_ptr<std::vector<std::string>> idxs(){ return idxs_;}
    std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
    std::shared_ptr<std::vector<bool>> aops(){ return aops_;}
    
    std::string Tsymm() {return Tsymm_; }
    std::string psymm(){ return psymm_; }
    
    std::shared_ptr<std::vector<int>> plus_ops(){ return plus_ops_;}
    std::shared_ptr<std::vector<int>> kill_ops(){ return kill_ops_;}
    
    std::vector<std::shared_ptr<TensOp<DataType>>> orig_tensors(){ return orig_tensors_ ;}
    std::shared_ptr<std::vector<int>> Tsizes() { return Tsizes_ ;}
    std::shared_ptr<std::vector<int>> cmlsizevec() { return cmlsizevec_ ; }
    std::shared_ptr<std::vector<std::string>> names() { return names_; }
    
    void generate_ranges();
    
    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors, std::shared_ptr<std::vector<std::string>> id_ranges );
    
    void get_ctrs_tens_ranges() override; 
    
    void print_gamma_contribs();

}; 
#endif
//     With spin construction 
//
//     MultiTensOp<DataType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
//                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
//                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};
