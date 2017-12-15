#ifndef __SRC_SMITH_TensOp_H
#define __SRC_SMITH_TensOp_H
#include <src/smith/wicktool/wickutils.h>
#include <src/smith/wicktool/ctrtensop.h>

 //#include "wickutils.h"
 //#include "ctrtensop.h"

//using namespace WickUtils;

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

class transformed_range_info {
    public :

    bool is_unique;
    std::shared_ptr<const std::vector<std::string>> unique_block;
    std::shared_ptr<const std::vector<std::string>> transformed_idxs;
    std::shared_ptr<const std::vector<std::pair<int,int>>> factors; 
                                
    transformed_range_info( bool is_unique_in,
                 std::shared_ptr<const std::vector<std::string>> unique_block_in,
                 std::shared_ptr<const std::vector<std::string>> transformed_idxs_in,
                 std::shared_ptr<const std::vector<std::pair<int,int>>> factors_in ) : 
                 is_unique(is_unique_in), unique_block(unique_block_in), transformed_idxs(transformed_idxs_in), factors(factors_in) {};
    ~transformed_range_info(){};

};


class TensOp_General {
   public:
     const std::vector<int> plus_ops_;
     std::shared_ptr<const std::vector<int>> plus_ops_ptr_;
  
     const std::vector<int> kill_ops_;
     std::shared_ptr<const std::vector<int>> kill_ops_ptr_;
    
     // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
     const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks_;
     std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks_ptr_;
  
     // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
     const std::map< const std::vector<std::string>,
                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >> all_ranges_;
 
     std::shared_ptr<const std::map< const std::vector<std::string>,
                                      std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>  >> all_ranges_ptr_;
  
  public:
    TensOp_General( std::vector<int> plus_ops, std::vector<int> kill_ops,
                    std::vector< std::shared_ptr<const std::vector<std::string>>> unique_range_blocks,
                    std::map< const std::vector<std::string>, 
                              std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>> all_ranges );
    ~TensOp_General(){};
  
    std::shared_ptr<const std::vector<int>> plus_ops() const { return plus_ops_ptr_; }
    std::shared_ptr<const std::vector<int>> kill_ops() const { return kill_ops_ptr_; }
    
    // unique_range_blocks tells you which parts of your tensor actually need to be calculated and stored.
    std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return  unique_range_blocks_ptr_;}
    
    std::shared_ptr<const std::map< const std::vector<std::string>,
                    std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                    all_ranges()const  {return  all_ranges_ptr_; }

};

class MultiTensOp_General : public  TensOp_General {
    public : 
    const std::vector<std::string> idxs_;
    const std::vector<bool> aops_;
    const std::vector<std::vector<std::string>> idx_ranges_;
    const std::pair<double,double> orig_factor_;
    const int num_idxs_;

    const std::vector<int> cmlsizevec_;

    std::shared_ptr<const std::vector<std::string>> idxs_ptr_;
    std::shared_ptr<const std::vector<bool>> aops_ptr_;
    std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges_ptr_;
    std::shared_ptr<const std::vector<int>> cmlsizevec_ptr_;

    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
    const std::map< const std::vector<std::string>,
    std::tuple<bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>>>> split_ranges_map_ ;

      MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                           std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                           std::vector< std::shared_ptr<const std::vector<std::string>>>& unique_range_blocks,
                           std::shared_ptr<std::map< const std::vector<std::string>, 
                           std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>> all_ranges, 
                           std::shared_ptr<std::map< const std::vector<std::string>,
                           std::tuple<bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >> split_ranges_map ) ; 
    ~MultiTensOp_General(){};

};


namespace TensOp {
template<typename DataType>
class TensOp {
   protected :  
     const std::string name_;

   private :
     const std::vector<std::string> idxs_;
     const std::vector<std::vector<std::string>> idx_ranges_;
     const std::vector<bool> aops_;
     const DataType orig_factor_;
     const std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
     const std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;
     const std::string Tsymm_;
     const int num_idxs_;
     bool apply_symmetry( const std::vector<std::string>& ranges_1, const std::vector<std::string>& ranges_2,
                          std::map< const std::vector<std::string>, 
                          std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int> > >& all_ranges_ );
     bool satisfies_constraints( std::vector<std::string>& ranges );
     std::shared_ptr<const TensOp_General> Op_dense_;

     std::shared_ptr< const std::vector<std::string>> idxs_const_ptr_; 
     std::shared_ptr< const std::vector<bool>> aops_const_ptr_; 
     std::shared_ptr< const std::vector<std::vector<std::string>>> idx_ranges_ptr_;
  

   public:
     TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
             std::vector<bool>& aops, DataType orig_factor,
             std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
             std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
             std::string& Tsymm );
     ~TensOp(){};

     int num_idxs(){ return num_idxs_; }

     int factor(){ return orig_factor_; };

     std::string const name(){ return name_;}

     std::string const Tsymm(){ return Tsymm_; }

     std::shared_ptr<const std::vector<std::string>> idxs(){ return idxs_const_ptr_;}

     std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_ptr_;}

     std::shared_ptr<const std::vector<bool>> aops(){ return aops_const_ptr_;}

     std::shared_ptr<const std::vector<int>> plus_ops(){ return Op_dense_->plus_ops();}

     std::shared_ptr<const std::vector<int>> kill_ops(){ return Op_dense_->kill_ops();}

     std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return Op_dense_->unique_range_blocks(); }
    
     std::shared_ptr<const std::map< const std::vector<std::string>,
                                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                                     all_ranges() const  {return Op_dense_->all_ranges(); }

     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     
     void generate_ranges();
     void get_ctrs_tens_ranges() ;

};
}


namespace MultiTensOp {
template<typename DataType>
class MultiTensOp  {

   private :
     std::shared_ptr<const MultiTensOp_General> Op_dense_;
     const std::string name_;

   public :
     std::vector<std::shared_ptr<TensOp::TensOp<DataType>>> orig_tensors_; 
    
     int num_tensors_;
     
     DataType orig_factor_;

     MultiTensOp( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>& orig_tensors );
    ~MultiTensOp(){};

    std::tuple<std::shared_ptr<std::map< const std::vector<std::string>, std::tuple<bool,                                      std::shared_ptr<const std::vector<std::string>>,                     std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >>,
               std::shared_ptr<std::map< const std::vector<std::string>, std::tuple<bool,                                      std::shared_ptr<const std::vector<std::string>>,                     std::shared_ptr<const std::vector<std::string>>, std::pair<int,int>>  >>, 
               std::shared_ptr<std::map< const std::vector<std::string>, std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >>>>
    generate_ranges( int num_idxs, std::vector<int>& cmlsizevec );
    void get_ctrs_tens_ranges(); 
    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors, std::shared_ptr<std::vector<std::string>> id_ranges );

    int num_idxs() { return Op_dense_->num_idxs_; }
    std::pair<double,double> factor() { return Op_dense_->orig_factor_; };

    std::string const  name()const { return name_;}
 
    std::shared_ptr<const std::vector<std::string>> idxs() const { return Op_dense_->idxs_ptr_; }

    std::shared_ptr<const std::vector<bool>> aops() const { return Op_dense_->aops_ptr_; }

    std::shared_ptr<const std::vector<int>> cmlsizevec() const { return Op_dense_->cmlsizevec_ptr_; }

    std::shared_ptr<const std::vector<int>> plus_ops() const { return Op_dense_->plus_ops_ptr_; }

    std::shared_ptr<const std::vector<int>> kill_ops() const { return Op_dense_->kill_ops_ptr_; }

    std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges() const { return Op_dense_->idx_ranges_ptr_; }

    std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return Op_dense_->unique_range_blocks(); }
    
    std::shared_ptr<const std::map< const std::vector<std::string>,
                                    std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                                    all_ranges() const  {return Op_dense_->all_ranges(); }

};
}
#endif
//     With spin construction 
//
//     MultiTensOp::MultiTensOp<DataType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
//                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
//                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};
