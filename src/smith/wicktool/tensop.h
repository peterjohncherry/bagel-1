#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
 #include <src/smith/wicktool/wickutils.h>
 #include <src/smith/wicktool/ctrtensop.h>

 //#include "wickutils.h"
 //#include "ctrtensop.h"

//using namespace WickUtils;

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

namespace TensOp_Interface {
template<typename DataType> 
 class TensOp_Interface {
    public:
      virtual ~TensOp_Interface(){};

    virtual int  num_idxs() = 0;
    virtual std::string  name() = 0;
    virtual std::string Tsymm() = 0;
    virtual DataType factor() = 0;
    virtual std::shared_ptr<const std::vector<std::string>> idxs() =0;
    virtual std::shared_ptr<const std::vector<bool>> aops() = 0;
    virtual std::shared_ptr<const std::vector<int>> plus_ops() = 0;
    virtual std::shared_ptr<const std::vector<int>> kill_ops() = 0;
    virtual std::shared_ptr<const std::map< const std::vector<std::string>,
                                            std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>> range_map() = 0;
    virtual std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() = 0 ;

};
}
namespace TensOp_General {
template<typename DataType> 
class TensOp_General {
   public:
     const std::vector<int> plus_ops_;
     const std::shared_ptr<const std::vector<int>> plus_ops_ptr_;

     const std::vector<int> kill_ops_;
     const std::shared_ptr<const std::vector<int>> kill_ops_ptr_;
    
     // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
     const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks_;
     const std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks_ptr_;

     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     const std::map< const std::vector<std::string>,
                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >> all_ranges_;
     const std::shared_ptr<const std::map< const std::vector<std::string>,
                                      std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>  >> all_ranges_ptr_;

  public:
    TensOp_General( std::vector<int> plus_ops, std::vector<int> kill_ops,
                    std::vector< std::shared_ptr<const std::vector<std::string>>> unique_range_blocks,
                    std::map< const std::vector<std::string>, 
                              std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>> all_ranges ) :
                    plus_ops_(plus_ops), plus_ops_ptr_(std::make_shared<const std::vector<int>>(plus_ops_)),
                    kill_ops_(kill_ops), kill_ops_ptr_(std::make_shared<const std::vector<int>>(kill_ops_)), 
                    unique_range_blocks_(unique_range_blocks), unique_range_blocks_ptr_(std::make_shared<const std::vector< std::shared_ptr< const std::vector<std::string>>>>(unique_range_blocks_)),
                    all_ranges_(all_ranges), all_ranges_ptr_( std::make_shared<const std::map< const std::vector<std::string>,
                                                                                              std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                                                              std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>(all_ranges_)) {};
    ~TensOp_General(){};

     std::shared_ptr<const std::vector<int>> plus_ops(){ return plus_ops_ptr_;}
     std::shared_ptr<const std::vector<int>> kill_ops(){ return kill_ops_ptr_;}
     
     // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
     std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks(){ return  unique_range_blocks_ptr_;}
     
     std::shared_ptr<const std::map< const std::vector<std::string>,
                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                     all_ranges() {return  all_ranges_ptr_; }

};
}

namespace TensOp {
template<typename DataType>
class TensOp {
   public:
     const std::string name_;
     const std::vector<std::string> idxs_;
     const std::vector<std::vector<std::string>> idx_ranges_;
     const std::vector<bool> aops_;
     const DataType orig_factor_;
     
     
     const std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
     const std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;
     const std::string Tsymm_;
     const int num_idxs_;

     //TensOp_General::TensOp_General<DataType> TensOp_Dense_Info;   
 
     std::shared_ptr< const std::vector<std::string>> idxs_const_ptr; 
        

     TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
                  std::vector<bool>& aops, DataType orig_factor,
                  std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
                  std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
                  std::string& Tsymm );
     ~TensOp(){};

     int num_idxs() { return 1 ; }
     int factor()   { return 1 ; };

     std::string const name(){ return  "a";}
     std::string const Tsymm(){ return "a"; }

     std::shared_ptr<const std::vector<std::string>> idxs(){ return std::make_shared<const std::vector<std::string>>( idxs_);}
     std::shared_ptr<const std::vector<bool>> aops(){ return std::make_shared<const std::vector<bool>>(aops_);}
//     std::shared_ptr<const std::vector<int>> plus_ops(){ return std::make_shared<const std::vector<int>>(plus_ops_);}
//     std::shared_ptr<const std::vector<int>> kill_ops(){ return std::make_shared<const std::vector<int>>(kill_ops_);}
     
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     
     void generate_ranges();
     bool apply_symmetry( const std::vector<std::string>& ranges_1, const std::vector<std::string>& ranges_2,
                          std::map< const std::vector<std::string>,  std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > >& all_ranges_ );
     bool satisfies_constraints( std::vector<std::string>& ranges );
     void get_ctrs_tens_ranges() ;

};
}


namespace MultiTensOp_Prep {
template<typename DataType>
class MultiTensOp_Prep {
   public:
     std::string name_;
     std::string Tsymm_;
     std::vector<std::string> idxs_;
     std::vector<std::vector<std::string>> idx_ranges_;
     std::vector<bool> aops_;
     
     std::vector<int> plus_ops_;
     std::vector<int> kill_ops_;
     int num_idxs_;

     int num_tensors_;
     
     DataType orig_factor_;

     std::vector<std::shared_ptr< const std::vector<std::string>>> split_idxs_ ;
     std::vector<std::string> names_  ;
     std::vector<int> cmlsizevec_ ;
     std::vector<int> Tsizes_ ;

     std::vector<std::shared_ptr<TensOp_Interface::TensOp_Interface<DataType>>> orig_tensors_; 

     // unique_range_blocks tells you what parts of your tensor actually need to be calculated and stored.
     std::vector< std::shared_ptr< const std::vector<std::string>> > unique_range_blocks;
     
     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     std::map< const std::vector<std::string>, 
               std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > > all_ranges_;

    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
    std::shared_ptr< std::map< const std::vector<std::string>, 
                               std::tuple<std::shared_ptr<const std::vector<bool>>, std::shared_ptr<const std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>>>>> split_ranges_map;

    //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
    std::shared_ptr<std::map< const std::vector<std::string>,
                              std::tuple<bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<const std::vector<std::pair<int,int>>> > >> combined_ranges;

    MultiTensOp_Prep( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp_Interface::TensOp_Interface<DataType>>>& orig_tensors );
    ~MultiTensOp_Prep(){};
     
    void generate_ranges();
    bool satisfies_constraints( std::vector<std::string>& ranges );
};
}

namespace MultiTensOp {
template<typename DataType>
//class MultiTensOp : public TensOp::TensOp<DataType> {
class MultiTensOp {
  public: 
      std::string name_;
      std::vector<std::string> idxs_;
      std::vector<std::vector<std::string>> idx_ranges_;
      std::vector<bool> aops_;
      DataType orig_factor_;
     
     
      std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
      std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;
      std::string Tsymm_;
      int num_idxs_;

    int num_idxs() { return 1 ; }
     int factor()   { return 1 ; };

     std::string  name(){ return  "a";}
     std::string  Tsymm(){ return "a"; }

     std::shared_ptr< std::vector<std::string>> idxs(){ return std::make_shared< std::vector<std::string>>( idxs_);}
     std::shared_ptr< std::vector<bool>> aops(){ return std::make_shared< std::vector<bool>>(aops_);}


     std::vector<std::shared_ptr< const std::vector<std::string>>> split_idxs_ ;
     std::vector<std::string> names_  ;
     std::vector<int> cmlsizevec_ ;
     std::vector<int> Tsizes_ ;

     std::vector<std::shared_ptr<TensOp_Interface::TensOp_Interface<DataType>>> orig_tensors_; 

     // unique_range_blocks tells you what parts of your tensor actually need to be calculated and stored.
     std::vector< std::shared_ptr< const std::vector<std::string>> > unique_range_blocks;
     
     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     std::map< const std::vector<std::string>, 
               std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > > all_ranges_;

    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
    std::shared_ptr< std::map< const std::vector<std::string>, 
                               std::tuple<std::shared_ptr<const std::vector<bool>>, std::shared_ptr<const std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>>>>> split_ranges_map;

    //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
    std::shared_ptr<std::map< const std::vector<std::string>,
                              std::tuple<bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<const std::vector<std::pair<int,int>>> > >> combined_ranges;













    // should be a better way; could define both MultiTensOp and TensOp as derived from non-templated class, as nothing but factors depends on DataType

    //cout maps from state and spin-sector to list of possible Aops
//    std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > ,
//                              std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths_;

//    std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors_;
//    bool spinfree_;
  
//    int num_tensors_;

//    std::shared_ptr< std::map< std::vector<std::string>, std::pair<int,int> > > orb_ranges;

  private:
    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
//    std::shared_ptr< std::map< std::vector<std::string>,
//                               std::tuple<std::shared_ptr<std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>>>>> split_ranges;

  public: 

    //Multi Tensor analogue of all_ranges; same as normal tensor but with vector of factors. The ordering of these factors is the same as that of the operators in the tensor
//    std::shared_ptr<std::map< std::vector<std::string>,
//                              std::tuple<bool, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> combined_ranges;

    //takes from pair of unc _range_ (defines gamma) and ctr_indexes (defines A contrib) 
//    std::shared_ptr<std::map< std::tuple< std::vector<std::string>, pstr_vec, pstr_vec >, std::shared_ptr<std::vector<std::string>> >> CMTP_gamma_contribs;
//    std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
    std::shared_ptr< std::map< std::string, std::shared_ptr<CtrMultiTensorPart<DataType>> >> CMTP_map; 

//    MultiTensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType> >> orig_tensors) :  TensOp::TensOp<DataType>::TensOp( name, spinfree, orig_tensors) {};
    MultiTensOp(std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType> >> orig_tensors) {};// :  TensOp::TensOp<DataType>::TensOp( name, spinfree, orig_tensors) {};
    ~MultiTensOp() {};
    
    void generate_ranges();
    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors, std::shared_ptr<std::vector<std::string>> id_ranges );
    void get_ctrs_tens_ranges(); 
    void print_gamma_contribs();
}; 
}
#endif
//     With spin construction 
//
//     MultiTensOp::MultiTensOp<DataType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
//                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
//                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};
