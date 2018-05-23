#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_ALGOP_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_ALGOP_H

//Classes for defining contraction operations
//ctr_abs pos is the position of the contracted index in the totally uncontracted tensor
//ctr_rel_pos is the position of the contracted index in the contracted tensor
#include <string> 
#include <vector> 
#include <utility> 
#include <memory>
#include <map>

class CtrOp_base {
  public : 
    const std::string Tout_name_;
    const std::string ctr_type_;

    CtrOp_base() {};
    CtrOp_base(std::string Tout_name_in, std::string ctr_type_in): Tout_name_(Tout_name_in), ctr_type_(ctr_type_in) {};
    ~CtrOp_base(){};

    virtual std::string Tout_name(){ return Tout_name_ ;};
    virtual std::string ctr_type(){ return ctr_type_ ;}
    virtual std::string T1name(){ throw std::runtime_error("Do not access T1name from CtrOp_base class!"); return "" ;}
    virtual std::string T2name(){ throw std::runtime_error("Do not access T2name from CtrOp_base class!"); return "";};
    virtual int T1_ctr_abs_pos(){ throw std::runtime_error("Do not access T1_ctr_abs_pos from CtrOp_base class!"); return 1;}; 
    virtual int T2_ctr_abs_pos(){ throw std::runtime_error("Do not access T2_ctr_abs_pos from CtrOp_base class!"); return 1;};
    virtual int T1_ctr_rel_pos(){ throw std::runtime_error("Do not access T1_ctr_rel_pos from CtrOp_base class!"); return 1;}; 
    virtual int T2_ctr_rel_pos(){ throw std::runtime_error("Do not access T2_ctr_rel_pos from CtrOp_base class!"); return 1;};
    virtual std::pair<int,int> ctr_rel_pos(){ throw std::runtime_error("Do not access ctr_rel_pos from CtrOp_base class!"); return std::make_pair(-1,-1) ;};
    virtual std::pair<int,int> ctr_abs_pos(){ throw std::runtime_error("Do not access ctr_abs_pos from CtrOp_base class!"); return std::make_pair(-1,-1) ;};;
    virtual std::shared_ptr<std::vector<int>> new_order() { throw std::runtime_error( "Do not access new_order from CtrOp_base class!");  return std::make_shared<std::vector<int>>(1,0);}
    
    virtual std::shared_ptr<std::vector<std::string>> tensor_list(){ throw std::runtime_error("Do not access tensor_list() from CtrOp_base class!");
              return std::make_shared<std::vector<std::string>>();
            };

};

 
class CtrOp_diff_T : public CtrOp_base {
  public : 
    const std::string T1name_;
    const std::string T2name_;
    const int T1_ctr_abs_pos_; 
    const int T2_ctr_abs_pos_;
    const int T1_ctr_rel_pos_; 
    const int T2_ctr_rel_pos_;

    CtrOp_diff_T(std::string T1name_in, std::string T2name_in , std::string Tout_name_in , int T1_ctr_abs_pos_in, int T2_ctr_abs_pos_in,
                 int T1_ctr_rel_pos_in, int T2_ctr_rel_pos_in, std::string ctr_type_in ): CtrOp_base(Tout_name_in, ctr_type_in),
    T1name_(T1name_in), T2name_(T2name_in) ,  T1_ctr_abs_pos_(T1_ctr_abs_pos_in),  T2_ctr_abs_pos_(T2_ctr_abs_pos_in), 
    T1_ctr_rel_pos_(T1_ctr_rel_pos_in),  T2_ctr_rel_pos_(T2_ctr_rel_pos_in) {};

    ~CtrOp_diff_T(){};
    std::string T1name() override { return T1name_ ;}
    std::string T2name() override { return T2name_ ;};
    int T1_ctr_abs_pos() override { return T1_ctr_abs_pos_;}; 
    int T2_ctr_abs_pos() override { return T2_ctr_abs_pos_;};
    int T1_ctr_rel_pos() override { return T1_ctr_rel_pos_;}; 
    int T2_ctr_rel_pos() override { return T2_ctr_rel_pos_;};

};
  
class CtrOp_same_T : public CtrOp_base {
  public : 
    const std::string T1name_;
    const std::pair<int,int> ctr_abs_pos_;
    const std::pair<int,int> ctr_rel_pos_;

    CtrOp_same_T( std::string T1name_in, std::string Tout_name_in, std::pair<int,int> ctr_abs_pos_in,
                  std::pair<int,int> ctr_rel_pos_in, std::string ctr_type_in ): CtrOp_base(Tout_name_in, ctr_type_in),
    T1name_(T1name_in), ctr_abs_pos_(ctr_abs_pos_in),  ctr_rel_pos_(ctr_rel_pos_in) {};

    ~CtrOp_same_T(){};

    std::string T1name() override { return T1name_ ;}
    std::pair<int,int> ctr_rel_pos() override { return ctr_rel_pos_;};
    std::pair<int,int> ctr_abs_pos() override { return ctr_abs_pos_;};
};      
 
 
class CtrOp_single_id : public CtrOp_base {
  public : 
    const std::string T1name_;
    const int T1_ctr_abs_pos_; 
    const int T1_ctr_rel_pos_; 

    CtrOp_single_id( std::string T1name, std::string T2name , std::string Tout_name , int T1_ctr_abs_pos, int T2_ctr_abs_pos,
                     int T1_ctr_rel_pos, int T2_ctr_rel_pos, std::string ctr_type ) : CtrOp_base(Tout_name, ctr_type),
                     T1name_(T1name), T1_ctr_abs_pos_(T1_ctr_abs_pos), T1_ctr_rel_pos_(T1_ctr_rel_pos) {};

    ~CtrOp_single_id(){};

    std::string T1name() override { return T1name_ ;}
    int T1_ctr_abs_pos() override { return T1_ctr_abs_pos_;}; 
    int T1_ctr_rel_pos() override { return T1_ctr_rel_pos_;}; 

};
 
class CtrOp_exc_ids : public CtrOp_base {
  public : 
    const std::string T1name_;
    const std::pair<int,int> ctr_abs_pos_;
    const std::pair<int,int> ctr_rel_pos_;


    CtrOp_exc_ids( std::string T1name, std::string Tout_name, std::pair<int,int> ctr_abs_pos,
                   std::pair<int,int> ctr_rel_pos, std::string ctr_type ) : CtrOp_base(Tout_name, ctr_type),
                   T1name_(T1name), ctr_abs_pos_(ctr_abs_pos),  ctr_rel_pos_(ctr_rel_pos) {};

    ~CtrOp_exc_ids(){};

    std::string T1name() override { return T1name_ ;}
    std::pair<int,int> ctr_rel_pos() override { return ctr_rel_pos_;};
    std::pair<int,int> ctr_abs_pos() override { return ctr_abs_pos_;};

};

//Just a reordering operation
//TODO find some way of combining this with the reordering from the gamma_map 
class CtrOp_reorder : public CtrOp_base {
  public : 
    const std::string T1name_;
    const std::pair<int,int> ctr_abs_pos_;
    const std::pair<int,int> ctr_rel_pos_;
    const std::shared_ptr<std::vector<int>> new_order_;

    CtrOp_reorder( std::string T1name, std::string Tout_name, std::shared_ptr<std::vector<int>> new_order, std::string ctr_type  ) :
                   CtrOp_base(Tout_name, ctr_type), T1name_(T1name), new_order_(new_order) {};

    ~CtrOp_reorder(){};

    std::string T1name() override { return T1name_ ;}
    std::shared_ptr<std::vector<int>> new_order() { return new_order_;}

};

// For performing direct products, but name is "cartesian" so as to avoid confusion with CtrOp_DiffT ...
class CtrOp_DirectProduct : public CtrOp_base {
  public : 
    std::shared_ptr<std::vector<std::string>> tensor_list_;

    CtrOp_DirectProduct( std::shared_ptr<std::vector<std::string>> tensor_list, std::string Tout_name, std::string ctr_type ) :
                          CtrOp_base(Tout_name, ctr_type), tensor_list_(tensor_list) {};

    ~CtrOp_DirectProduct(){};

    std::shared_ptr<std::vector<std::string>> tensor_list() { return tensor_list_; }

};

//Just tells it to get the tensor; better to have all operations, even ones this simple, spat out from the compute list
//TODO currently does not need to be seperate class to base, but for now keep seperate for consistencies sake 
class CtrOp_Get : public CtrOp_base {
  public : 

    CtrOp_Get( std::string Tout_name, std::string ctr_type ) : CtrOp_base(Tout_name, ctr_type) {};

    ~CtrOp_Get(){};

};
 
#endif
