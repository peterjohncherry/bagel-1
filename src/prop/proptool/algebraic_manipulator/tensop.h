#ifndef __SRC_PROP_PROPTOOL_TensOp_H
#define __SRC_PROP_PROPTOOL_TensOp_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

class Op_General_base { 

     protected:
       const std::vector<std::string> idxs_;
       const std::vector<bool> aops_;
       const std::vector<int> plus_ops_;
       const std::vector<int> kill_ops_;
       const std::vector<std::vector<std::string>> idx_ranges_;
       const std::pair<double,double> orig_factor_;
       const int num_idxs_;

     public:
       const std::vector< std::shared_ptr< const std::vector<std::string>>> unique_range_blocks_;
     
       std::shared_ptr<const std::vector<std::string>> idxs_ptr_;
       std::shared_ptr<const std::vector<bool>> aops_ptr_;
       
       std::shared_ptr<const std::vector<int>> plus_ops_ptr_;
       std::shared_ptr<const std::vector<int>> kill_ops_ptr_;
       
       std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges_ptr_;
       std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks_ptr_;

       std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<range_block_info>>> all_ranges_;
       std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<split_range_block_info>>> split_ranges_;

      public :

        Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                        std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                        std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks, 
                        std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<range_block_info>>> all_ranges);


        Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                        std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                        std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks, 
                        std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<split_range_block_info>>> split_ranges );
       ~Op_General_base(){};
       
        int num_idxs() const { return num_idxs_; }
       
	std::pair<double,double> factor() const { return orig_factor_; };

        std::shared_ptr<const std::vector<std::string>> idxs() const { return idxs_ptr_;}
        std::string idxs(int ii ) const { return idxs_ptr_->at(ii);}
       
        std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges() const { return idx_ranges_ptr_;}
        std::vector<std::string>  idx_ranges(int ii ) const { return idx_ranges_ptr_->at(ii);}
       
        std::shared_ptr<const std::vector<bool>> aops() const { return aops_ptr_;}
        bool aops(int ii) const { return aops_ptr_->at(ii) ;}
       
        std::shared_ptr<const std::vector<int>> plus_ops()const { return plus_ops_ptr_;}
        int plus_ops(int ii ) const { return plus_ops_ptr_->at(ii) ;}
       
        std::shared_ptr<const std::vector<int>> kill_ops() const{ return kill_ops_ptr_;}
        int kill_ops(int ii ) const { return kill_ops_ptr_->at(ii) ;}
       
        std::shared_ptr<const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return  unique_range_blocks_ptr_;}
        std::shared_ptr< const std::vector<std::string>> unique_range_blocks(int ii) const { return  unique_range_blocks_ptr_->at(ii);}
       
      // add in more virtual functions for range blocks and state dependence

        std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<range_block_info>>> all_ranges() const  {return  all_ranges_; }
        std::shared_ptr<range_block_info> all_ranges(const std::vector<std::string> range_block ) const  { return all_ranges_->at(range_block) ; }

        virtual std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info>>> split_ranges()const { return split_ranges_; }
        virtual std::shared_ptr<split_range_block_info> split_ranges(const std::vector<std::string> range_block )const { return split_ranges_->at(range_block) ; }

        virtual std::shared_ptr<const std::vector<int>> cmlsizevec() const  = 0 ;
        virtual int cmlsizevec(int ii) const  = 0 ;

};

class TensOp_General : public Op_General_base {

  public:

    TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                    std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                    std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                    std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr< range_block_info> >> all_ranges ); 

    TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                    std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                    std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                    std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info >>> all_ranges);  
 
   ~TensOp_General(){};

    std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info > > > split_ranges() const{
      throw std::logic_error("1 Should not be trying to access split range map from merged TensOp_General, probably error in range looping! Aborting! " ) ; 
    return split_ranges_; } 

    std::shared_ptr< split_range_block_info > split_ranges(const std::vector<std::string> range_block )const {
      throw std::logic_error("2 Should not be trying to access split range map from merged TensOp_General, probably error in range looping! Aborting! " ) ; 
    return split_ranges_->at(range_block); } 

    std::shared_ptr<const std::vector<int>> cmlsizevec() const  { std::cout << "WARNING!! this is a single tensor.." << std::endl; return std::make_shared<std::vector<int>>(0,0); } ;
    int cmlsizevec(int ii )const { if ( ii > 0 ) throw std::logic_error( "this is a single tensor...., something wrong in TensOp") ; return 0;}

};

class MultiTensOp_General : public  Op_General_base {

    private:
    
    const std::vector<int> cmlsizevec_;
    std::shared_ptr<const std::vector<int>> cmlsizevec_ptr_;

    public : 
      MultiTensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                           std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor, std::vector<int>& cmlsizevec, 
                           std::shared_ptr<std::vector< std::shared_ptr<const std::vector<std::string>>>> unique_range_blocks,
                           std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<split_range_block_info>>> all_ranges );
    ~MultiTensOp_General(){};

     std::shared_ptr<const std::vector<int>> cmlsizevec() const { return cmlsizevec_ptr_;}
     int cmlsizevec(int ii )const { return cmlsizevec_[ii];}

};


class TensOp_base {

   protected :
     const std::string name_;
     const bool spinfree_;
     std::string Tsymm_;
     int state_dep_; 
     std::shared_ptr<const Op_General_base> Op_dense_;

   public:

 
     TensOp_base( std::string name, bool spinfree, std::string Tsymm, int state_dep ) : name_(name), spinfree_(spinfree), Tsymm_(Tsymm), state_dep_(state_dep)  {};
     TensOp_base( std::string name, bool spinfree ) : name_(name), spinfree_(spinfree), Tsymm_("none"), state_dep_(0)  {};
     ~TensOp_base(){};

     std::string const name(){ return name_;}

     std::string const Tsymm(){ return Tsymm_;}
     
     int num_idxs(){ return Op_dense_->num_idxs();}

     std::pair<double,double> factor(){ return Op_dense_->factor(); };

     std::shared_ptr< const std::vector<std::string>> idxs(){ return Op_dense_->idxs();}
     std::string idxs(int ii ){ return Op_dense_->idxs(ii);}

     std::shared_ptr< const std::vector<bool>> aops(){ return Op_dense_->aops();}
     bool aops(int ii){ return Op_dense_->aops(ii);}

     std::shared_ptr< const std::vector<std::vector<std::string>>> idx_ranges(){ return Op_dense_->idx_ranges();}
     std::vector<std::string> idx_ranges(int ii){ return Op_dense_->idx_ranges(ii);}

     std::shared_ptr< const std::vector<int>> plus_ops(){ return Op_dense_->plus_ops();}
     int plus_ops(int ii){ return Op_dense_->plus_ops(ii);}

     std::shared_ptr< const std::vector<int>> kill_ops(){ return Op_dense_->kill_ops();}
     int kill_ops(int ii){ return Op_dense_->kill_ops(ii);}


     std::shared_ptr< const std::vector< std::shared_ptr< const std::vector<std::string>>>> unique_range_blocks() const { return Op_dense_->unique_range_blocks(); }
     std::shared_ptr< const std::vector<std::string>> unique_range_blocks(int ii) const { return Op_dense_->unique_range_blocks(ii); }
     
     virtual std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<range_block_info > > > all_ranges() const  { return Op_dense_->all_ranges(); }
     virtual std::shared_ptr< range_block_info > all_ranges(const std::vector<std::string> range_block ) const  { return Op_dense_->all_ranges(range_block); }

     virtual std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info > > > split_ranges() const{ return Op_dense_->split_ranges(); } 
     virtual std::shared_ptr< split_range_block_info > split_ranges(const std::vector<std::string> range_block )const { return Op_dense_->split_ranges(range_block); } 
//     virtual std::shared_ptr<range_block_info > sparsed_ranges(const std::vector<int> range_block ) const  = 0 ; 

};

namespace TensOp {
template<typename DataType>
class TensOp : public TensOp_base {
   private :
     std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > > symmfuncs_; 
     std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) > constraints_;

     bool satisfies_constraints( std::vector<std::string>& ranges );

//   protected : 
//     std::map< std::vector<int>, std::vector<std::shared_ptr<const std::vector<std::string> > >  > state_sparsity_map_; 

   public:
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart_Base> > > CTP_map_;

     TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
             std::vector<bool>& aops, DataType orig_factor,
             std::vector< std::tuple< std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>), int, int > >& symmfuncs, 
             std::vector<bool(*)(std::shared_ptr<std::vector<std::string>>) >& constraints,
             std::string& Tsymm, int state_dep);
     ~TensOp(){};

     void get_ctrs_tens_ranges() ;

     void get_ctp_idxs_ranges( std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos, std::shared_ptr<range_block_info> block_info );

     std::tuple< std::shared_ptr< const std::map< const std::vector<std::string> , std::shared_ptr<range_block_info> >>,
                                  std::shared_ptr<std::vector< std::shared_ptr< const std::vector<std::string> >>> >
     generate_ranges( std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges, std::vector<bool>& aops );
  
     std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<range_block_info > > > all_ranges() const  { return Op_dense_->all_ranges(); }
     std::shared_ptr< range_block_info > all_ranges(const std::vector<std::string> range_block ) const  { return Op_dense_->all_ranges(range_block); }

     std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info > > > split_ranges() const{
       throw std::logic_error("1 TensOp Should not be trying to access split range map from merged TensOp, probably error in range looping! Aborting! " ) ; 
     return Op_dense_->split_ranges(); } 

     std::shared_ptr< split_range_block_info > split_ranges(const std::vector<std::string> range_block )const {
       throw std::logic_error("2 TensOp Should not be trying to access split range map from merged TensOp, probably error in range looping! Aborting! " ) ; 
     return Op_dense_->split_ranges(range_block); } 

    std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base > >> CTP_map() { return CTP_map_; } 
};
}

namespace MultiTensOp {
template<typename DataType>
class MultiTensOp : public TensOp_base {

   public :
     std::vector<std::shared_ptr<TensOp::TensOp<DataType>>> orig_tensors_; 

     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart_Base> > > CTP_map_;

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> > >> CMTP_map_;

     int num_tensors_;
     
     DataType orig_factor_;

     MultiTensOp( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp::TensOp<DataType>>>& orig_tensors );
    ~MultiTensOp(){};

    std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info >>>
    generate_ranges( std::vector<std::string>& idxs, std::vector<bool>& aops, std::vector<int>& cmlsizevec );

    std::shared_ptr<const std::vector<int>> cmlsizevec() const  {  return Op_dense_->cmlsizevec(); } ;
    int cmlsizevec(int ii )const { return Op_dense_->cmlsizevec(ii); };
 
    void get_ctrs_tens_ranges(); 
  
    std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<range_block_info > > > all_ranges() const  {
      return std::dynamic_pointer_cast<const std::map< const std::vector<std::string>, std::shared_ptr<range_block_info > > >(Op_dense_->all_ranges());
    }

    std::shared_ptr< range_block_info > all_ranges(const std::vector<std::string> range_block ) const  { 
      return std::dynamic_pointer_cast<range_block_info>(Op_dense_->split_ranges(range_block));
    }

    std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<split_range_block_info > > > split_ranges() const{return Op_dense_->split_ranges(); } 
    std::shared_ptr< split_range_block_info > split_ranges(const std::vector<std::string> range_block )const { return Op_dense_->split_ranges(range_block); } 

    void enter_into_CMTP_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges );

    std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base> >> CTP_map() { return CTP_map_; } 
    std::shared_ptr< std::map< std::string, std::shared_ptr< CtrMultiTensorPart<DataType> > >> CMTP_map() { return CMTP_map_; } 
};
}
#endif
