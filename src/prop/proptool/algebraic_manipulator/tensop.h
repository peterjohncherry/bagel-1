#ifndef __SRC_PROP_PROPTOOL_TensOp_H
#define __SRC_PROP_PROPTOOL_TensOp_H
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/ctrtensop.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/symmetry_operations.h>
#include <src/prop/proptool/algebraic_manipulator/constraints.h>

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
     std::shared_ptr<std::set<std::shared_ptr<Range_Block_Info>>> required_blocks_;

     // MultiTens specific, but should be obtainable for TensOp (with warning)
     std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_; 


   public:

     std::pair<double,double> factor_;

     std::shared_ptr<std::map<std::string, std::shared_ptr<std::map<const std::vector<std::string>, std::shared_ptr<Range_Block_Info>>>>> all_ranges_state_specific_;
    
     TensOp_Base( std::string name, bool spinfree, std::string Tsymm, int state_dep ) : name_(name),  factor_(std::make_pair(1.0,1.0)), spinfree_(spinfree),
                                                                                        Tsymm_(Tsymm), state_dep_(state_dep),
                                                                                        required_blocks_(std::make_shared<std::set<std::shared_ptr<Range_Block_Info>>>()) {;
                                                                                        std::cout << " TensOp_Base cstr1" << std::endl; };

     TensOp_Base( std::string name, std::pair<double,double>& factor, bool spinfree ) : name_(name), factor_(factor), spinfree_(spinfree), Tsymm_("none"), state_dep_(0),
                                                                                        required_blocks_(std::make_shared<std::set<std::shared_ptr<Range_Block_Info>>>()){
                                                                                        std::cout << " TensOp_Base cstr2" << std::endl; };

     TensOp_Base( std::string name, bool spinfree, std::vector<std::shared_ptr<TensOp_Base>>& sub_tensops );

     ~TensOp_Base(){};

     std::string name() const { return name_;}

     std::string Tsymm() const { return Tsymm_;}
     
     int num_idxs() const { return num_idxs_;}

     std::pair<double,double> factor() const { return factor_; };

     void print_definition() const;

     std::shared_ptr< const std::vector<std::string>> idxs() const { return idxs_;}
     std::string idxs(int ii ){ return (*idxs_)[ii];}

     std::shared_ptr< const std::vector<bool>> aops() const { return aops_;}
     bool aops(int ii){ return (*aops_)[ii];}

     std::shared_ptr< const std::vector<std::vector<std::string>>> idx_ranges(){ return idx_ranges_;}
     std::vector<std::string> idx_ranges(int ii){ return (*idx_ranges_)[ii];}

     virtual bool satisfies_constraints( std::vector<std::string>& ranges ) { return true ; } 
    
     void add_required_block( std::shared_ptr<Range_Block_Info> block ) { required_blocks_->emplace( block ); } 

     std::shared_ptr< std::map< std::string, std::shared_ptr< CtrTensorPart_Base> >> CTP_map() { return CTP_map_; } 

     virtual void generate_ranges( std::shared_ptr<Op_Info> op_info ) = 0;

     virtual void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); }

     virtual void add_state_ids( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); } 

     virtual std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> state_ids() {
                throw std::logic_error("should not call from anything but TensOp");
                std::shared_ptr<std::set<std::shared_ptr<Op_Info>>> dummy;
                return dummy; } 

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const char op_trans ) {
         throw std::logic_error("Should not call from base class"); return std::make_shared<std::vector<bool>>(0); }

     virtual std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ) = 0;

     virtual bool is_projector(){ return false ; } 

     virtual std::shared_ptr< std::map< const std::vector<std::string>, std::shared_ptr<Range_Block_Info > > > all_ranges() const  { return all_ranges_; }

     virtual void generate_uncontracted_ctps( std::shared_ptr<Op_Info> op_info ) {} ;

     virtual void generate_uncontracted_ctps( std::vector<std::vector<int>>& state_ids ) {} ;

     virtual void enter_cmtps_into_map(const pint_vec& ctr_pos_list, const std::vector<std::string>& id_ranges, std::shared_ptr<Op_Info> op_info ) = 0 ; 

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

    void apply_symmetry( const std::vector<std::string>& new_block, const std::vector<bool>& aops, std::shared_ptr<Op_Info> op_info );

    void add_to_range_block_map( std::vector<std::string>& idx_ranges );
 
    std::shared_ptr< std::map< char, std::tuple< std::pair<double,double>, std::shared_ptr<std::vector<std::vector<std::string>>> >>>   idx_ranges_map_;

  public:

   TensOp( std::string name, std::vector<std::string>& idxs, std::vector<std::vector<std::string>>& idx_ranges,
           std::vector<bool>& aops, std::pair<double,double>& orig_factor,
           std::vector<std::shared_ptr<Transformation>>& symmfuncs,
           std::vector<std::shared_ptr<Constraint>>& constraints,
           std::string& Tsymm, int state_dep, std::shared_ptr<std::map<char, long unsigned int>> range_prime_map  );
   ~TensOp(){};

   void generate_uncontracted_ctps( std::shared_ptr<Op_Info> op_info );

   void generate_blocks();

   void generate_ranges( std::shared_ptr<Op_Info> op_info );

   void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info );

   std::shared_ptr<Range_Block_Info> 
   get_transformed_range_block( std::shared_ptr<Op_Info> op_info, std::shared_ptr<Range_Block_Info>& block );

   void enter_cmtps_into_map(const pint_vec& ctr_pos_list, const std::vector<std::string>& id_ranges, std::shared_ptr<Op_Info> op_info ); 

   std::vector<std::shared_ptr<TensOp_Base>> sub_tensops(){  std::vector<std::shared_ptr<TensOp_Base>> sub_tensops_ ; return sub_tensops_; } 
  
   std::shared_ptr<std::vector<bool>> transform_aops( const char trans ); 

   std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

   bool satisfies_constraints( std::vector<std::string>& ranges ); 

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
 
    void generate_uncontracted_ctps( std::shared_ptr<Op_Info> op_info );

    void get_cmtp( std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>  ctp_vec, 
                   std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>>> ccp_vec,
                   std::shared_ptr<Op_Info> op_info );

    void shift_ccp_and_ctp_vecs( std::shared_ptr<CtrMultiTensorPart<DataType>>& tatb_cmtp,
                                 int ta, int tb, std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>>& ctp_vec,
                                 std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >>& ccp_vec  );

    void enter_cmtps_into_map(const pint_vec& ctr_pos_list, const std::vector<std::string>& id_ranges, std::shared_ptr<Op_Info> op_info );

    std::shared_ptr<std::vector<bool>> transform_aops( const std::vector<int>& op_order , const std::vector<char>& op_trans ); 

    virtual void generate_transformed_ranges( std::shared_ptr<Op_Info> op_info ) { throw std::logic_error("should not call from anything but TensOp"); }

};
}
#endif
