#ifndef __SRC_PROP_PROPTOOL_TensOp_H
#define __SRC_PROP_PROPTOOL_TensOp_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>


// TODO : FIX THIS STUPID STRUCTURE; I SUGGEST THE FOLLOWING:
//DQ : class A {//stuff      } 
//     class B : A { // more stuff }   
//     class C { A name }
//     class D : C {  B name }
//     c.f. Op_Dense. Would prefer to this without having function in D. 

//     Following solution looks OK:

//     class X<param>;
//     template specialization 
//     
//     class X<C> = A, class  X<D> = B;
//     class Y<C> = C, class  Y<D> = D;
//
//     class Y<C> { X<C> name } 
//     class Y<D> { X<D> name } 

// Instead of having a parent and load of virtuals, you have a base class, and specializations for the relevant templates.
// Is this way better? Presumably one has an advantage...
// It doesn't look it to me; you still end up with a four classes, and then a templated map, which you probably want to factor out.... so six again :(
// However, whilst you do technically have the same number of classes, you can at least have it so 

class TensOp_Base;
namespace TensOp {  template<typename DataType> class TensOp; } 
namespace MultiTensOp { template<typename DataType> class MultiTensOp; } 

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

class Op_General_base { 
     friend TensOp_Base;
     friend TensOp::TensOp<double>;
     friend TensOp::TensOp<std::complex<double>>;
     friend MultiTensOp::MultiTensOp<double>;
     friend MultiTensOp::MultiTensOp<std::complex<double>>;
   
     protected:
       const std::vector<std::string> idxs_;
       const std::vector<bool> aops_;
       const std::vector<int> plus_ops_;
       const std::vector<int> kill_ops_;
       const std::vector<std::vector<std::string>> idx_ranges_;
       const std::pair<double,double> orig_factor_;
       const int num_idxs_;
     
       const bool hermitian_ =  true; //TODO have this set, although should default to true.
       const bool hermitian_inverse_ =  true; //TODO have this set, although should default to false; true for projectors.

     public:
     
       std::shared_ptr<const std::vector<std::string>> idxs_ptr_;
       std::shared_ptr<const std::vector<bool>> aops_ptr_;
       
       std::shared_ptr<const std::vector<int>> plus_ops_ptr_;
       std::shared_ptr<const std::vector<int>> kill_ops_ptr_;
       
       std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges_ptr_;

       std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>> all_rxnges_;

       std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_rxnges_;
 
      public :

        Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                         std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                         std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<Range_Block_Info>>> all_rxnges);

        Op_General_base( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                         std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                         std::shared_ptr<const std::map< const std::vector<std::string>,  std::shared_ptr<Split_Range_Block_Info>>> split_rxnges );

       ~Op_General_base(){};
       
        int num_idxs() const { return num_idxs_; }
       
	std::pair<double,double> factor() const { return orig_factor_; };

        std::shared_ptr<const std::vector<std::string>> idxs() const { return idxs_ptr_;}
        std::string idxs(int ii ) const { return idxs_[ii];}
       
        std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges() const { return idx_ranges_ptr_;}
        std::vector<std::string>  idx_ranges(int ii ) const { return idx_ranges_[ii];}
       
        std::shared_ptr<const std::vector<bool>> aops() const { return aops_ptr_;}
        bool aops(int ii) const { return aops_[ii] ;}
       
        std::shared_ptr<const std::vector<int>> plus_ops()const { return plus_ops_ptr_;}
        int plus_ops(int ii ) const { return plus_ops_[ii] ;}
       
        std::shared_ptr<const std::vector<int>> kill_ops() const{ return kill_ops_ptr_;}
        int kill_ops(int ii ) const { return kill_ops_[ii] ;}

        bool hermitian() { return hermitian_; }  

        // add in more virtual functions for range blocks and state dependence
        std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>> all_rxnges() const  {return  all_rxnges_; }
    

        virtual std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_rxnges()const { return split_rxnges_; }

        virtual std::shared_ptr<const std::vector<int>> cmlsizevec() const  = 0 ;
        virtual int cmlsizevec(int ii) const  = 0 ;

};

class TensOp_General : public Op_General_base {

  public:

    TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                    std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                    std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr< Range_Block_Info > >> all_rxnges ) :   
                    Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, all_rxnges) {};

    TensOp_General( std::vector<std::string>& idxs,  std::vector<bool>& aops, std::vector<int>& plus_ops, std::vector<int>& kill_ops,
                    std::vector<std::vector<std::string>>& idx_ranges, std::pair<double,double> factor,
                    std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info >>> split_rxnges ) : 
                    Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, split_rxnges ) {};
 
   ~TensOp_General(){};

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
                           std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr< Split_Range_Block_Info>>> split_rxnges ) :
                           Op_General_base::Op_General_base( idxs, aops, plus_ops, kill_ops, idx_ranges, factor, split_rxnges), 
                           cmlsizevec_(cmlsizevec), cmlsizevec_ptr_(std::make_shared<const std::vector<int>>(cmlsizevec_)) {};

     ~MultiTensOp_General(){};

      std::shared_ptr<const std::vector<int>> cmlsizevec() const { return cmlsizevec_ptr_;}
      int cmlsizevec(int ii )const { return cmlsizevec_[ii];}

};


class TensOp_Base {


   protected :

     const std::string name_;
     const bool spinfree_;
     std::string Tsymm_;
     int state_dep_; 
     std::shared_ptr<const Op_General_base> Op_dense_;
     std::shared_ptr<std::set<std::string>> required_blocks_;
     std::shared_ptr< std::map< std::string, std::shared_ptr<CtrTensorPart_Base> > > CTP_map_;
     std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_; 
     bool projector_; 

     std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>> hconj_ranges_;


   public:
 
     TensOp_Base( std::string name, bool spinfree, std::string Tsymm, int state_dep ) : name_(name), spinfree_(spinfree), Tsymm_(Tsymm), state_dep_(state_dep),
                                                                                        required_blocks_(std::make_shared<std::set<std::string>>()) {};

     TensOp_Base( std::string name, bool spinfree ) : name_(name), spinfree_(spinfree), Tsymm_("none"), state_dep_(0), 
                                                      required_blocks_(std::make_shared<std::set<std::string>>())  {};

     TensOp_Base( std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops );

     ~TensOp_Base(){};

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

     virtual 
     bool satisfies_constraints( std::vector<std::string>& ranges ) { return true ; } 
    
     void add_required_block( std::string block_name ) { required_blocks_->emplace( block_name ); } 
     std::shared_ptr<std::set<std::string>> required_blocks( std::string block_name ) { return required_blocks_; } 

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base> >> CTP_map() { return CTP_map_; } 

     void transform_aops_rngs( std::vector<char>& rngs, std::pair<double,double>& factor, const char op_trans_in ); 
  
     void generate_transformed_ranges( char transformation );

     virtual
     std::shared_ptr<std::vector<char>>
     transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::pair<double, double>& factor, 
                          const std::vector<int>& op_order , const std::vector<char>& op_trans ) {
       throw  std::logic_error( " should not be in split transform aops_rngs in base class " );
       std::shared_ptr<std::vector<char>> dummy;
       return dummy;
     } 
           
     std::shared_ptr<Range_Block_Info> transform_block_rngs( const std::vector<char>& rngs, const char op_trans_in );

     virtual  
     std::shared_ptr<Split_Range_Block_Info> 
     transform_block_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::shared_ptr<std::vector<bool>> trans_aops,
                           const std::vector<int>& op_order, const std::vector<char>& op_trans ){
       throw  std::logic_error( " should not be in split transform aops_rngs in base class " );
       return block;
     } 

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const char op_trans ) = 0;

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ) = 0;

     virtual bool is_projector(){ return false ; } 

     virtual std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info > > > all_rxnges() const  { return Op_dense_->all_rxnges(); }

     virtual void generate_uncontracted_ctps() = 0;

     virtual std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info > > > split_rxnges() const = 0 ;  

     virtual void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges ) = 0 ; 
     virtual std::vector<std::shared_ptr<TensOp_Base>> sub_tensops() = 0;
 
};

namespace TensOp {
template<typename DataType>
class TensOp :  public TensOp_Base , public std::enable_shared_from_this<TensOp<DataType>> {

   private :
     
     std::vector<std::function<void( std::vector<std::string>& )>> symmfuncs_;
     std::vector<std::function<bool( std::vector<std::string>& )>> constraints_;

     void add_to_range_block_map( std::vector<std::string>& idx_ranges );

     std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >  all_ranges_tmp_;

     std::shared_ptr<const std::vector<std::string>> 
     apply_symmetry( std::vector<int>& idxs_trans, const std::vector<std::string>& new_block, std::pair<double,double>& fac_new );


 
   public:

     TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
             std::vector<bool>& aops, DataType orig_factor,
             std::vector<std::function<void( std::vector<std::string>& )>>& symmfuncs_,
             std::vector<std::function<bool( std::vector<std::string>& )>>& constraints_,
             std::string& Tsymm, int state_dep, std::shared_ptr<std::map<char, long unsigned int>> range_prime_map  );
     ~TensOp(){};

     void generate_uncontracted_ctps();

     std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>
     generate_rangesX( std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges, std::vector<bool>& aops );

    std::shared_ptr<const std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_rxnges() const{return Op_dense_->split_rxnges(); } 

     void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges ) { 
        throw std::logic_error( "TensOp::TensOp<DataType> should cannot call enter_into_CTP_map form this class!! Aborting!!" ); } 

    std::vector<std::shared_ptr<TensOp_Base>> sub_tensops(){  std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_ ; return sub_tensops_; } 
   
    std::shared_ptr<std::vector<bool>> transform_aops( const char trans ); 

    std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

    void transform_aops_rngs( std::vector<char>& rngs, std::pair<double,double>& factor, const char op_trans_in ); 
  
    std::shared_ptr<std::vector<char>>
    transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, const std::vector<int>& op_order , const std::vector<char>& op_trans );

    bool satisfies_constraints( std::vector<std::string>& ranges ); 
};
}

namespace MultiTensOp {
template<typename DataType>
class MultiTensOp : public TensOp_Base, public std::enable_shared_from_this<MultiTensOp<DataType>> {

   public :

     int num_tensors_;
     
     DataType orig_factor_;

     MultiTensOp( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops,
                  std::shared_ptr<std::map< char , long unsigned int >> range_prime_map );

    ~MultiTensOp(){};

    std::shared_ptr< const std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info> >>
    generate_rangesX( std::vector<std::string>& idxs, std::vector<bool>& aops, std::vector<int>& cmlsizevec  ); 

    std::shared_ptr<const std::vector<int>> cmlsizevec() const  {  return Op_dense_->cmlsizevec(); } ;
    int cmlsizevec(int ii )const { return Op_dense_->cmlsizevec(ii); };

    std::vector<std::shared_ptr<TensOp_Base>> sub_tensops(){ return sub_tensops_; } 
 
    void generate_uncontracted_ctps();

    void get_cmtp( std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>  ctp_vec, 
                   std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>> ccp_vec );

    void shift_ccp_and_ctp_vecs( std::shared_ptr<CtrMultiTensorPart<DataType>>& tatb_cmtp,
                                 int ta, int tb, std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>& ctp_vec,
                                 std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >>& ccp_vec  );

    std::shared_ptr<const std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_rxnges() const{return Op_dense_->split_rxnges(); } 

    void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges );

    std::shared_ptr<std::vector<bool>>  transform_aops( const char trans ); 

    std::shared_ptr<std::vector<bool>>  transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

    std::shared_ptr<std::vector<char>> transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::pair<double, double>& factors,
                                                            const std::vector<int>& op_order, const std::vector<char>& op_trans );

    std::shared_ptr<Split_Range_Block_Info> 
    transform_block_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::shared_ptr<std::vector<bool>> trans_aops, 
                          const std::vector<int>& op_order, const std::vector<char>& op_trans ); 

};
}
#endif
