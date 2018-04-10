#ifndef __SRC_PROP_PROPTOOL_WICKTOOL_TERM_INIT_H
#define __SRC_PROP_PROPTOOL_WICKTOOL_TERM_INIT_H

#include <src/prop/proptool/proputils.h>
//Classes for initializing objects in equation
class Op_Init {
  
  private :
    const std::vector<std::string> idxs_;
    const std::shared_ptr<std::vector<int*>> idx_ptrs_;
  
  public :
    std::string name_;
    std::string alg_name_;
    int state_dep_ = 0;
    bool proj_op_;
    char trans_; // H or h : hermitian conjugate, T or t = transpose, C or c = complex conjugate, anything else : nothing                     
     
    Op_Init(std::string base_name, std::vector<std::string>& idxs, std::shared_ptr<std::vector<int*>> idx_ptrs ) :
            name_(base_name), idxs_(idxs), idx_ptrs_(idx_ptrs), alg_name_(base_name), trans_('n') {
              if (idxs.front() != "none" ) {
                alg_name_ += "_{";
                for (std::string idx : idxs_ )
                  alg_name_ += idx;
                alg_name_ += "}";
                state_dep_ = idxs_.size();
              } else {
                state_dep_ = 0;
              }
            };
 
    Op_Init( std::string base_name, std::vector<std::string>& idxs, std::shared_ptr<std::vector<int*>> idx_ptrs,
             std::string trans ) :
             name_(base_name), idxs_(idxs), idx_ptrs_(idx_ptrs), alg_name_(base_name), trans_(trans[0]) {
              std::cout << "trans = " << trans << std::endl;
              std::cout << "trans_ = " << trans_ << std::endl;
              if (idxs.front() != "none" ) {
                alg_name_ += "_{";
                for (std::string idx : idxs_ )
                  alg_name_ += idx;
                alg_name_ += "}";
                state_dep_ = idxs_.size();
              } else {
                state_dep_ = 0;
              }
            };

   ~Op_Init(){};

     int num_idxs(){ return idxs_.size(); }

     std::string state_specific_name(){
       std::string ssname = name_ ;
       ssname = + "_{";
       for ( auto id_val : *idx_ptrs_ )
         ssname += std::to_string(*id_val);
       ssname +="}" ;
       return ssname;
     }

     std::shared_ptr<std::vector<int>> get_idxs() {
       std::shared_ptr<std::vector<int>> state_idxs_list = std::make_shared<std::vector<int>>(idxs_.size());
       if (idxs_.front() != "none" ) {
         for ( int ii = 0 ; ii != state_idxs_list->size(); ii++ )
         state_idxs_list->at(ii) = *(idx_ptrs_->at(ii));
       }
       return state_idxs_list;
     }

     void get_op_idxs(std::vector<int>& state_idxs_list ) { // to avoid incredible number of shared_ptrs
       if (idxs_.front() != "none" ) {
         for ( int ii = 0 ; ii != state_idxs_list.size(); ii++ )
           state_idxs_list[ii] = *(idx_ptrs_->at(ii));
       }
       return;
     }

};

// BraKet e.g. < M | HT^{LN} | N >
// note, at present the code cannot handle summations in BraKets, e.g.,
// < M | H (f-E_{L}+E_{S})T^{LN} | N > must be specified as
// < M | HfT^{LN} | N > + E_{L}< M | HT^{LN} | N > +E_{s}< M | HT^{LN} | N >
// though these terms will still be merged (if possible) in expression
class BraKet_Init {

   std::string bra_index_;
   int* bra_index_ptr_;

   std::string ket_index_;
   int* ket_index_ptr_;

  public :
   std::string name_;
   std::shared_ptr<std::vector<Op_Init>> op_list_;
   std::string proj_name_;
   bool projected_;

   BraKet_Init( std::shared_ptr<std::vector<Op_Init>> op_list,
                std::string bra_index, int* bra_index_ptr, std::string ket_index, int* ket_index_ptr ) :
                op_list_(op_list), bra_index_(bra_index), bra_index_ptr_(bra_index_ptr),
                ket_index_(ket_index),  ket_index_ptr_(ket_index_ptr), projected_(false) {
               
                  name_ = "<" + bra_index + "|";
                  for ( Op_Init op : *op_list )
                    name_ += op.alg_name_;
                  
                  name_+= "|" + ket_index_+ ">" ;
               };

   BraKet_Init( std::shared_ptr<std::vector<Op_Init>> op_list,
                std::string bra_index, int* bra_index_ptr, std::string ket_index, int* ket_index_ptr, std::string proj_name ) :
                op_list_(op_list), bra_index_(bra_index), bra_index_ptr_(bra_index_ptr),
                ket_index_(ket_index),  ket_index_ptr_(ket_index_ptr), proj_name_(proj_name), projected_(true) {
               
                  name_ = "<" + bra_index + "|";
                  for ( Op_Init op : *op_list )
                    name_ += op.alg_name_;
                  
                  name_+= "|" + ket_index_+ ">" ;
               };




   ~BraKet_Init(){};

   std::shared_ptr<std::vector<std::vector<int>>> get_op_idxs_list() {
     std::shared_ptr<std::vector<std::vector<int>>> op_idxs_list = std::make_shared<std::vector<std::vector<int>>>(op_list_->size());
     for ( int ii = 0 ; ii != op_idxs_list->size(); ii++ ){
       op_idxs_list->at(ii) = std::vector<int>( op_list_->at(ii).num_idxs() );
       op_list_->at(ii).get_op_idxs( op_idxs_list->at(ii) );
     }
     return op_idxs_list;
   }

   void get_op_idxs_list( std::vector<std::vector<int>>& op_idxs_list ) {
     for ( int ii = 0 ; ii != op_idxs_list.size(); ii++ ){
       op_idxs_list[ii] = std::vector<int>( op_list_->at(ii).num_idxs() );
       op_list_->at(ii).get_op_idxs( op_idxs_list[ii] );
     }
     return;
   }
    
//   Op_Init& op_list( int jj ) { return op_list_->at(jj); }
   int bra_index(){ return *bra_index_ptr_; }
   int ket_index(){ return *ket_index_ptr_; }
   

};
//Just a list of BraKets, this class appears superflous, but is necessary for handling of summations.
//All indices will be summed  over simultaneously, e.g.,
// can do  :
//      sum_{M} <M| H | N>
// but not  :
//   \sum_{L} ( E_{L}+sum_{M} <M| H | N> )
// as the latter could need to sum over N and _then_ L .
class Term_Init {

  public :
    std::string name_;
    std::string type_;
    std::shared_ptr<std::vector<BraKet_Init>> braket_list_;
    std::shared_ptr<std::vector<std::string>> braket_factors_;
    std::shared_ptr<std::map<std::string, int>> idx_val_map_;
    std::shared_ptr<std::map<std::string, std::string>> idx_name_map_;
    bool orbital_projector_;
    std::string proj_op_name_;

    std::string alg_name_;

    Term_Init( std::string name, std::string type,
               std::shared_ptr<std::vector<BraKet_Init>> braket_list,
               std::shared_ptr<std::vector<std::string>> braket_factors,
               std::shared_ptr<std::map<std::string, int>> idx_val_map) :
               name_(name), type_(type), braket_list_(braket_list), braket_factors_(braket_factors),
               idx_val_map_(idx_val_map), orbital_projector_(false) {
             
               alg_name_ = "";
               for ( int ii =0 ; ii != braket_factors_->size(); ii++ )
                 alg_name_ += "(" + braket_factors_->at(ii) + ")" + braket_list->at(ii).name_ + " + ";
               alg_name_.pop_back();
               std::cout << "======================= New Term =======================" << std::endl;
               std::cout << alg_name_ << std::endl << std::endl;
               };
    
    // second variation of initialization for terms with orbital projector
    Term_Init( std::string name, std::string type,
               std::shared_ptr<std::vector<BraKet_Init>> braket_list,
               std::shared_ptr<std::vector<std::string>> braket_factors,
               std::shared_ptr<std::map<std::string, int>> idx_val_map,
               std::string proj_op_name ) :
               name_(name), type_(type), braket_list_(braket_list), braket_factors_(braket_factors),
               idx_val_map_(idx_val_map), orbital_projector_(true), proj_op_name_(proj_op_name)  {
             
               alg_name_ = "";
               for ( int ii =0 ; ii != braket_factors_->size(); ii++ )
                 alg_name_ += "(" + braket_factors_->at(ii) + ")" + braket_list->at(ii).name_ + " + ";
               alg_name_.pop_back();
               std::cout << "======================= New Term =======================" << std::endl;
               std::cout << alg_name_ << std::endl << std::endl;
               };

    ~Term_Init(){};
};
// A list of terms, all terms (i.e. all BraKets ) will be merged
// Different terms can have different terms, but note that the index variables are globally defined for the expression.
// Care must be taken so that indices are not inadvertently tied to together, e.g.,
// If we want to evaluate
// X = sum_{M} <M| H | N> +  sum_{L}<M| E_{L} | N>
// And I define Expression  using terms
// Term 1 = sum_{M} <M| H |N>
// Term 2 = sum_{L} <M| E_{L} | Q>
// Then the program will actually evaluate
// Term 1 = sum_{M} <M| H | N> +  sum_{M}sum_{L}<M| E_{L} | N>
// i.e., it will sum over M in the second term, which is incorrect.
// To avoid this one can just specify different index names, e.g.,
// Term 1 = sum_{M} <M| H |N>
// Term 2 = sum_{L} <R| E_{L} |N>
// Note that once the expression has been initialized all variable names vanish and are replaced by indexes, hence
// this need to change index names should not (I think) impact braket/term/expression reuse.
class Expression_Init {

  public :
    std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init>>>> term_list_;

    //key :: index name
    //result :: pair ( true_if_we_sum_over_this_index ,  name_of_index_range ) 
    std::shared_ptr<std::vector<std::shared_ptr<std::map<std::string,std::pair<bool,std::string>>>>> term_range_maps_;
   
    std::string name_;
    
    std::string type_;

    Expression_Init( std::shared_ptr<std::vector<std::pair<std::string, std::shared_ptr<Term_Init>>>> term_list,
                     std::shared_ptr<std::vector<std::shared_ptr<std::map<std::string,std::pair<bool,std::string>>>>> term_range_maps,
                     std::string type ):
                     term_list_(term_list), term_range_maps_(term_range_maps), type_(type) {
                       for ( std::pair<std::string, std::shared_ptr<Term_Init>> term : *term_list_ )
                         name_ += "(" +term.first +")."+ term.second->name_+ "+";
                       name_.pop_back();
                       std::cout << "Expression name " << std::endl;
                       std::cout << name_ << std::endl << std::endl;
                     };
    ~Expression_Init(){};
 
};
#endif
