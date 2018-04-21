#ifndef __SRC_PROP_PROPTOOL_TensOp_H
#define __SRC_PROP_PROPTOOL_TensOp_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>

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

class TensOp_Base {

   protected :
    
     // fields which should be set on initialization
     const std::string name_;
     std::shared_ptr<const std::vector<std::string>> idxs_;
     std::shared_ptr<const std::vector<bool>> aops_;
     std::shared_ptr<const std::vector<std::vector<std::string>>> idx_ranges_;
     std::pair<double,double> factor_;

     //TODO keep like this for now
     const bool projector_ = false;
     const bool spinfree_  = false;
     const bool hermitian_ =  true;
     const bool hermitian_inverse_ =  true;
     std::string Tsymm_;
     int state_dep_;

     // fields determined in constructor, which should not change during lifetime of class instance
     int num_idxs_;

     std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>> all_ranges_;

     std::shared_ptr<std::map<std::string, std::shared_ptr<CtrTensorPart_Base>>> CTP_map_;
 
     // TensOp specfic, but should be obtainable for MultiTensOp (with warning)
     std::shared_ptr<const std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>> hconj_ranges_;
     std::shared_ptr<std::set<std::string>> required_blocks_;

     // MultiTens specific, but should be obtainable for TensOp (with warning)
     std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_; 
     std::shared_ptr< std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_ranges_;

   public:

     std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>>>> all_ranges_state_specific_;
    
     std::shared_ptr<std::map<std::string,std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>>>> state_specific_split_ranges_; 
 
     TensOp_Base( std::string name, bool spinfree, std::string Tsymm, int state_dep ) : name_(name), spinfree_(spinfree), Tsymm_(Tsymm), state_dep_(state_dep),
                                                                                        required_blocks_(std::make_shared<std::set<std::string>>()) {};

     TensOp_Base( std::string name, std::pair<double,double>& factor, bool spinfree ) : name_(name), factor_(factor), spinfree_(spinfree), Tsymm_("none"),
                                                                                        state_dep_(0), required_blocks_(std::make_shared<std::set<std::string>>()) {};

     TensOp_Base( std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops );

     ~TensOp_Base(){};

     std::string const name(){ return name_;}

     std::string const Tsymm(){ return Tsymm_;}
     
     int num_idxs(){ return num_idxs_;}

     std::pair<double,double> factor(){ return factor_; };

     std::shared_ptr< const std::vector<std::string>> idxs(){ return idxs_;}
     std::string idxs(int ii ){ return (*idxs_)[ii];}

     std::shared_ptr< const std::vector<bool>> aops(){ return aops_;}
     bool aops(int ii){ return (*aops_)[ii];}

     std::shared_ptr< const std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
     std::vector<std::string> idx_ranges(int ii){ return (*idx_ranges_)[ii];}

     virtual 
     bool satisfies_constraints( std::vector<std::string>& ranges ) { return true ; } 
    
     void add_required_block( std::string block_name ) { required_blocks_->emplace( block_name ); } 
     std::shared_ptr<std::set<std::string>> required_blocks( std::string block_name ) { return required_blocks_; } 

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base> >> CTP_map() { return CTP_map_; } 

     void transform_aops_rngs( std::vector<char>& rngs, std::pair<double,double>& factor, const char op_trans_in ); 
  
     virtual void generate_ranges( std::shared_ptr<Op_Info> op_info ) = 0;

     virtual void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); }

     virtual void add_state_ids( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); } 

     virtual std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> state_ids() {
                throw std::logic_error("should not call from anything but TensOp");
                std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> dummy;
                return dummy; } 

     virtual
     std::shared_ptr<std::vector<char>>
     transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::pair<double, double>& factor, 
                          const std::vector<int>& op_order , const std::vector<char>& op_trans ) {
       throw  std::logic_error( " should not be in split transform aops_rngs in base class " );
       std::shared_ptr<std::vector<char>> dummy;
       return dummy;
     } 

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const char op_trans ) = 0;

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ) = 0;

     virtual bool is_projector(){ return false ; } 

     virtual std::shared_ptr< std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info > > > all_ranges() const  { return all_ranges_; }

     virtual std::shared_ptr< std::map< const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info > > > split_ranges() const = 0 ;  

     virtual void generate_uncontracted_ctps( std::shared_ptr<Op_Info> op_info ) {} ;

     virtual void generate_uncontracted_ctps( std::vector<std::vector<int>>& state_ids ) {} ;

     virtual void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges ) = 0 ; 

     virtual std::vector<std::shared_ptr<TensOp_Base>> sub_tensops() = 0;
 
};

namespace TensOp {
template<typename DataType>
class TensOp : public TensOp_Base , public std::enable_shared_from_this<TensOp<DataType>> {

  private :
    
    std::vector<std::shared_ptr<Transformation>> symmfuncs_;

    std::vector<std::shared_ptr<Constraint>> constraints_;

    std::vector<std::shared_ptr<std::vector<std::string>>> block_list_;
  
    std::shared_ptr<std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info> >>  all_ranges_tmp_;

    std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> state_ids_;

    void apply_symmetry( const std::vector<std::string>& new_block, const std::vector<bool>& aops, std::shared_ptr<Op_Info> op_info );

    void add_to_range_block_map( std::vector<std::string>& idx_ranges );
 
    std::shared_ptr< std::map< char, std::tuple< std::pair<double,double>, std::shared_ptr<std::vector<std::vector<std::string>>> >>>   idx_ranges_map_;

    void generate_idx_ranges_map();

  public:

   TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
           std::vector<bool>& aops, std::pair<double,double>& orig_factor,
           std::vector<std::shared_ptr<Transformation>>& symmfuncs,
           std::vector<std::shared_ptr<Constraint>>& constraints,
           std::string& Tsymm, int state_dep, std::shared_ptr<std::map<char, long unsigned int>> range_prime_map  );
   ~TensOp(){};

   void generate_uncontracted_ctps( std::shared_ptr<Op_Info> state_ids );

   void generate_blocks();

   void generate_ranges( std::shared_ptr<Op_Info> op_info );

   void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info );

   std::shared_ptr<Range_Block_Info> 
   get_transformed_range_block( std::shared_ptr<Op_Info> op_info, std::shared_ptr<Range_Block_Info>& block );

   std::shared_ptr< std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_ranges() const{ return split_ranges_; } 

   void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges ) { 
     throw std::logic_error( "TensOp::TensOp<DataType> should cannot call enter_into_CTP_map form this class!! Aborting!!" ); } 

   std::vector<std::shared_ptr<TensOp_Base>> sub_tensops(){  std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_ ; return sub_tensops_; } 
  
   std::shared_ptr<std::vector<bool>> transform_aops( const char trans ); 

   std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

   void transform_aops_rngs( std::vector<char>& rngs, std::pair<double,double>& factor, const char op_trans_in ); 
  
   std::shared_ptr<std::vector<char>>
   transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, const std::vector<int>& op_order , const std::vector<char>& op_trans );

   bool satisfies_constraints( std::vector<std::string>& ranges ); 

   void add_state_ids( std::shared_ptr<Op_Info> op_info ) { state_ids_->emplace(op_info); return; }

   std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> state_ids() { return state_ids_; } 

   std::shared_ptr<const std::vector<std::string>>
   apply_direct_range_transformation( const std::vector<std::string>& block, std::pair<double,double>& factor,
                                      std::shared_ptr<Op_Info>& op_info  );

};
}

namespace MultiTensOp {
template<typename DataType>
class MultiTensOp : public TensOp_Base, public std::enable_shared_from_this<MultiTensOp<DataType>> {

   public :

     int num_tensors_;

     std::shared_ptr<const std::vector<int>> cmlsizevec_;

     MultiTensOp( std::string name , bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops,
                  std::shared_ptr<std::map< char , long unsigned int >> range_prime_map );

    ~MultiTensOp(){};

    void generate_ranges( std::shared_ptr<Op_Info> op_info ); 

    std::shared_ptr<const std::vector<int>> cmlsizevec() const  {  return cmlsizevec_; } ;

    std::vector<std::shared_ptr<TensOp_Base>> sub_tensops(){ return sub_tensops_; } 
 
    void generate_uncontracted_ctps( std::shared_ptr<MultiOp_Info> op_info );

    void get_cmtp( std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>  ctp_vec, 
                   std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>> ccp_vec );

    void shift_ccp_and_ctp_vecs( std::shared_ptr<CtrMultiTensorPart<DataType>>& tatb_cmtp,
                                 int ta, int tb, std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>& ctp_vec,
                                 std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >>& ccp_vec  );

    std::shared_ptr< std::map<const std::vector<std::string>, std::shared_ptr<Split_Range_Block_Info>>> split_ranges() const { return split_ranges_; } 

    void enter_cmtps_into_map(pint_vec ctr_pos_list, std::pair<int,int> ReIm_factors, const std::vector<std::string>& id_ranges );

    std::shared_ptr<std::vector<bool>> transform_aops( const char trans ); 

    std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

    std::shared_ptr<std::vector<char>> transform_aops_rngs( std::shared_ptr<Split_Range_Block_Info> block, std::pair<double, double>& factors,
                                                            const std::vector<int>& op_order, const std::vector<char>& op_trans );

    virtual void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); }

};
}
#endif
