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

    private :
    const bool is_unique_;
    const bool survives_;
    const std::pair<int,int> factors_; 
    const std::shared_ptr<const std::vector<std::string>> orig_block_;
    const std::shared_ptr<const std::vector<std::string>> unique_block_;
    const std::shared_ptr<const std::vector<std::string>> transformed_idxs_;
                                
    public :
    transformed_range_info( bool is_unique,
                            bool survives,
                            std::pair<int,int> factors, 
                            std::shared_ptr<const std::vector<std::string>> orig_block,
                            std::shared_ptr<const std::vector<std::string>> unique_block,
                            std::shared_ptr<const std::vector<std::string>> transformed_idxs) :
                            orig_block_(orig_block), is_unique_(is_unique),
                            survives_(survives), unique_block_(unique_block), transformed_idxs_(transformed_idxs), factors_(factors) {};
    ~transformed_range_info(){};

    bool is_unique() const { return is_unique_ ; }
    bool survives() const { return survives_ ; }
    std::pair<int,int> factors() const { return factors_; }  
    std::shared_ptr<const std::vector<std::string>> orig_block()       { return orig_block_;      }
    std::shared_ptr<const std::vector<std::string>> unique_block()     { return unique_block_;    }
    std::shared_ptr<const std::vector<std::string>> transformed_idxs() { return transformed_idxs_;}

};


class TensOp_General {
   private:

     const std::vector<std::string> idxs_;
     const std::vector<bool> aops_;
     const std::vector<int> plus_ops_;
     const std::vector<int> kill_ops_;
     const std::vector<std::vector<std::string>> idx_ranges_;
     const std::pair<double,double> orig_factor_;
     const int num_idxs_;
     const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks_;
     const std::map< const std::vector<std::string>,
                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >> all_ranges_;
     const std::map< const std::vector<std::string>, std::shared_ptr<const transformed_range_info> > all_ranges_new_; 

   public:
     std::shared_ptr<const std::vector<std::string>> idxs_ptr_;
     std::shared_ptr<const std::vector<bool>> aops_ptr_;
     
     std::shared_ptr<const std::vector<int>> plus_ops_ptr_;
     std::shared_ptr<const std::vector<int>> kill_ops_ptr_;
     
     std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges_ptr_;
     std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks_ptr_;
     std::shared_ptr<const std::map< const std::vector<std::string>,
                                       std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,
                                                   std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>  >> all_ranges_ptr_;
  
  public:
    TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                    std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                    std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                    std::shared_ptr<std::map< const std::vector<std::string>, 
                    std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>> all_ranges); 
    ~TensOp_General(){};

     int num_idxs() const { return num_idxs_; }

     std::pair<double,double> factor() const { return orig_factor_; };

     std::shared_ptr<const std::vector<std::string>> idxs() const { return idxs_ptr_;}

     std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges() const { return idx_ranges_ptr_;}

     std::shared_ptr<const std::vector<bool>> aops() const { return aops_ptr_;}

     std::shared_ptr<const std::vector<int>> plus_ops()const { return plus_ops_ptr_;}

     std::shared_ptr<const std::vector<int>> kill_ops() const{ return kill_ops_ptr_;}

    std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return  unique_range_blocks_ptr_;}
    
    // all_ranges takes a possible rangeblock, and maps it to a unique rangeblock(1), a list of indexes(2)  and a factor(3)  resulting from the symmetry transformation
    std::shared_ptr<const std::map< const std::vector<std::string>,
                    std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                    all_ranges()const  {return  all_ranges_ptr_; }

};

class MultiTensOp_General : public  TensOp_General {
    public : 

    const std::vector<int> cmlsizevec_;

    //Similar to all_ranges, but has the ranges associated with the different operators split up, this may be useful for decomposition, so keep for now
    const std::map< const std::vector<std::string>,
    std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >> split_ranges_map_;

    std::shared_ptr<const std::vector<int>> cmlsizevec_ptr_;

    std::shared_ptr<std::map< const std::vector<std::string>,
    std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >>> split_ranges_map_ptr_;

    MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                         std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                         std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                         std::shared_ptr<std::map< const std::vector<std::string>, 
                         std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>> all_ranges, 
                         std::shared_ptr<std::map< const std::vector<std::string>,
                         std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >>> split_ranges_map );
    ~MultiTensOp_General(){};

     std::shared_ptr<const std::vector<int>> cmlsizevec() const { return cmlsizevec_ptr_;}
     int cmlsizevec(int ii )const { return cmlsizevec_[ii];}

//     std::shared_ptr<std::map< const std::vector<std::string>,
//     std::tuple<bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >>  split_ranges_map() const { return split_ranges_map_ptr_;}

};


namespace TensOp {
template<typename DataType>
class TensOp {
   protected :  
     const std::string name_;
     const bool spinfree_;

   private :
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;
     std::string Tsymm_;

     bool satisfies_constraints( std::vector<std::string>& ranges );

     std::shared_ptr<const TensOp_General> Op_dense_;

   public:
     TensOp( std::string name, bool spinfree ) : name_(name), spinfree_(spinfree){};
     TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
             std::vector<bool>& aops, DataType orig_factor,
             std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
             std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
             std::string& Tsymm );
     ~TensOp(){};

     std::string const name(){ return name_;}

     std::string const Tsymm(){ return Tsymm_;}

     int num_idxs(){ return Op_dense_->num_idxs();}

     std::pair<double,double> factor(){ return Op_dense_->factor(); };

     std::shared_ptr<const std::vector<std::string>> idxs(){ return Op_dense_->idxs();}

     std::shared_ptr<const std::vector<bool>> aops(){ return Op_dense_->aops();}

     std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges(){ return Op_dense_->idx_ranges();}

     std::shared_ptr<const std::vector<int>> plus_ops(){ return Op_dense_->plus_ops();}

     std::shared_ptr<const std::vector<int>> kill_ops(){ return Op_dense_->kill_ops();}

     std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return Op_dense_->unique_range_blocks(); }
    
     std::shared_ptr<const std::map< const std::vector<std::string>,
                                     std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                                     all_ranges() const  {return Op_dense_->all_ranges(); }

     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> > > CTP_map ;
     

     std::tuple< std::shared_ptr<std::map< const std::vector<std::string>, 
            std::tuple< bool, std::shared_ptr<const std::vector<std::string>>, std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>>>>,
            std::shared_ptr<std::vector< std::shared_ptr< const std::vector<std::string>> >>> generate_ranges( int num_idxs, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges );
     void get_ctrs_tens_ranges() ;

};
}

namespace MultiTensOp {
template<typename DataType>
class MultiTensOp : public  TensOp::TensOp<DataType> {

   private :
     const std::string name_;
     const bool spinfree_ = true;
     std::shared_ptr<const MultiTensOp_General> Op_dense_;

   public :
     std::vector<std::shared_ptr<TensOp::TensOp<DataType>>> orig_tensors_; 
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart<DataType>> >> CTP_map;
     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> > >> CMTP_map;
     int num_tensors_;
     
     DataType orig_factor_;

     MultiTensOp( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>& orig_tensors );
    ~MultiTensOp(){};

    std::tuple<std::shared_ptr<std::map< const std::vector<std::string>, std::tuple<bool,                                      std::shared_ptr<const std::vector<std::string>>,                     std::shared_ptr<const std::vector<std::string>>, std::shared_ptr<std::vector<std::pair<int,int>>> > >>,
               std::shared_ptr<std::map< const std::vector<std::string>, std::tuple<bool,                                      std::shared_ptr<const std::vector<std::string>>,                     std::shared_ptr<const std::vector<std::string>>, std::pair<int,int>>  >>, 
               std::shared_ptr<std::map< const std::vector<std::string>, std::tuple< std::shared_ptr<const std::vector<bool>>, std::shared_ptr<std::vector<std::shared_ptr<const std::vector<std::string>>>>, std::shared_ptr<std::vector<std::pair<int,int>>> >>>>
    generate_ranges( int num_idxs, std::vector<int>& cmlsizevec );

    std::shared_ptr<const std::vector<int>> cmlsizevec() const { return Op_dense_->cmlsizevec();}

    int cmlsizevec(int ii ) const { return Op_dense_->cmlsizevec(ii) ;}

    std::string const name(){ return name_;}

    int num_idxs(){ return Op_dense_->num_idxs();}

    std::pair<double,double> factor(){ return Op_dense_->factor(); };

    std::shared_ptr<const std::vector<std::string>> idxs(){ return Op_dense_->idxs();}

    std::shared_ptr<const std::vector<bool>> aops(){ return Op_dense_->aops();}

    std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges(){ return Op_dense_->idx_ranges();}

    std::shared_ptr<const std::vector<int>> plus_ops(){ return Op_dense_->plus_ops();}

    std::shared_ptr<const std::vector<int>> kill_ops(){ return Op_dense_->kill_ops();}

    std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return Op_dense_->unique_range_blocks(); }
    
    std::shared_ptr<const std::map< const std::vector<std::string>,
                                    std::tuple< bool, std::shared_ptr<const std::vector<std::string>>,  std::shared_ptr< const std::vector<std::string>>, std::pair<int,int>  >>>
                                    all_ranges() const  { return Op_dense_->all_ranges(); }

    void get_ctrs_tens_ranges(); 

    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges );

};
}
#endif
//     With spin construction 
//
//     MultiTensOp::MultiTensOp<DataType>(std::string name, std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors,
//                 std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths)
//                 : name_(name), spin_sectors_(spin_sectors), spin_paths_(spin_paths), spinfree_(false) {};
