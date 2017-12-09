#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/ctrtensop.h>

 //#include "wickutils.h"
 //#include "ctrtensop.h"

using namespace WickUtils;

template<class DataType>
class TensOp {
   protected:
   /* const */std::string Tsymm_;
   /* const */std::string psymm_;
   /* const */std::string name_;
   /* const */std::shared_ptr<std::vector<std::string>> idxs_;
   /* const */std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges_;
   /* const */std::shared_ptr<std::vector<bool>> aops_;
               
   /* const */std::shared_ptr<std::vector<int>> plus_ops_;
   /* const */std::shared_ptr<std::vector<int>> kill_ops_;
   /* const */std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    plus_combs_;
   /* const */std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    kill_combs_;

   private:
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_;
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints_;  
  
   /*  const*/ DataType factor_;

     std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;
     std::shared_ptr< std::map< std::vector<std::string>, 
                      std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges_;
     
     void get_ctrs_pos() ;
     //map from original index range, to ranges and factors used in calculation
     //tuple contained 1: bool is_this_range_unique?, 2: unique range which needs to be calculated, 3: indexes for contraction of this range,  4: factor from transformation 

   public:
     TensOp() {};
     TensOp( std::string name,
             std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs, 
             std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints) :
             name_(name), symmfuncs_(symmfuncs), constraints_(constraints), factor_(1.0) {};

     ~TensOp(){};
  
     std::string name(){ return name_;}
     std::shared_ptr<std::vector<std::string>> idxs(){ return idxs_;}
     std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
     std::shared_ptr<std::vector<bool>> aops(){ return aops_;}
     
     std::string Tsymm() {return Tsymm_; }
     std::string psymm(){ return psymm_; }


    std::shared_ptr<std::vector<int>> plus_ops(){ return plus_ops_;}
    std::shared_ptr<std::vector<int>> kill_ops(){ return kill_ops_;}
    std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > > plus_combs(){ return plus_combs_;}
    std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > > kill_combs(){ return kill_combs_;}

 
     int factor() { return factor_; };

     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> >,
                               std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >>spin_paths(){ return spin_paths_;}

     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs()
     { return symmfuncs_;}

     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints(){ return constraints_;}


     bool contracted;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
 
     virtual void initialize( std::vector<std::string> orig_idxs, std::vector<std::vector<std::string>> orig_idx_ranges,
                              std::vector<bool> orig_aops, int orig_factor, std::string orig_Tsymm, std::string orig_psymm = "2el" );
     
     virtual // change the mapped type so is [urpose built class; want to avoid tuples for the sake of clarity.
     std::shared_ptr<std::map< std::vector<std::string>, 
                     std::tuple<bool, std::shared_ptr<std::vector<std::string>>,  std::shared_ptr<std::vector<std::string>>, std::pair<int,int> > >> all_ranges(){ return all_ranges_; }

     //should be set up to generate hermitian conjugated operator TODO check to see if this is working
     void hconj();
  
     void generate_ranges( std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges, std::string symm_type);
     virtual void get_ctrs_tens_ranges() ;

     bool apply_symmetry( std::shared_ptr<std::vector<std::string>> ranges_1, std::shared_ptr<std::vector<std::string>> ranges_2  );
     bool apply_symmetry( std::vector<std::string> ranges_1, std::vector<std::string> ranges_2  );

     bool satisfies_constraints( std::shared_ptr<std::vector<std::string>> ranges );
     bool satisfies_constraints( std::vector<std::string> ranges );

     void generate_all_contractions();
 
}; 

template<class DataType>
class MultiTensOp : public TensOp<DataType> {
   public: 

     std::string name_;
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_;
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>)> constraints_;  
     std::shared_ptr<std::vector<std::string>> idxs;
     std::shared_ptr<std::vector<std::vector<std::string>>> idx_ranges;
     std::shared_ptr<std::vector<bool>> aops;
     std::shared_ptr<DataType> data;
     bool contracted;
     int factor;

     std::shared_ptr<std::vector<int>> plus_ops;
     std::shared_ptr<std::vector<int>> kill_ops;
     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    plus_combs;
     std::shared_ptr< std::vector< std::shared_ptr< std::vector< std::shared_ptr<std::vector<int> > > > > >    kill_combs;

     std::vector<std::shared_ptr<TensOp<DataType>>> orig_tensors_;
     std::shared_ptr<std::vector<int>> Tsizes ;
     std::shared_ptr<std::vector<int>> cmlsizevec ;
     std::shared_ptr<std::vector<std::string>> names;

     std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> split_idxs   ;
     std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> split_plus_ops;
     std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> split_kill_ops;

     //cout maps from state and spin-sector to list of possible Aops
     std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > ,
                               std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;

     std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors_;
     bool spinfree_;
  
     int num_tensors;

     std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrMultiTensorPart<DataType>> >> CMTP_map; 

     //takes from pair of unc _range_ (defines gamma) and ctr_indexes (defines A contrib) 
     std::shared_ptr<std::map< std::tuple< std::vector<std::string>, pstr_vec, pstr_vec >, std::shared_ptr<std::vector<std::string>> >> CMTP_gamma_contribs;

   private:
     //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
     std::shared_ptr< std::map< std::vector<std::vector<std::string>>, 
                                std::tuple<std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>>>>> split_ranges;

   public: 
     MultiTensOp() {};
     MultiTensOp(std::string name, bool spinfree) : name_(name), spinfree_(spinfree) {}; // spin free construction

     ~MultiTensOp() {};

     //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
     std::shared_ptr<std::map< std::vector<std::string>,
                               std::tuple<bool, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> combined_ranges;
    
    void initialize(std::vector<std::shared_ptr<TensOp<DataType> >> orig_tensors);

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
