#ifndef __SRC_SMITH_CtrTensOp_H
#define __SRC_SMITH_CtrTensOp_H

  #include <src/smith/wicktool/WickUtils.h>
//  #include "WickUtils.h"

//Classes for defining contraction operations
//ctr_abs pos is the position of the contracted index in the totally uncontracted tensor
//ctr_rel_pos is the position of the contracted index in the contracted tensor
class CtrOp_base {
  public : 
    const std::string Tout_name_;
    const std::string ctr_type_;

    CtrOp_base() {};
    CtrOp_base(std::string Tout_name_in, std::string ctr_type_in): Tout_name_(Tout_name_in), ctr_type_(ctr_type_in) {};
    ~CtrOp_base(){};

    virtual std::string Tout_name(){ return Tout_name_ ;};
    virtual std::string ctr_type(){ return ctr_type_ ;}
    virtual std::string T1name() { throw std::runtime_error("Not defined in CtrOp_base class!"); return "" ;}
    virtual std::string T2name(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return "";};
    virtual int T1_ctr_abs_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return 1;}; 
    virtual int T2_ctr_abs_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return 1;};
    virtual int T1_ctr_rel_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return 1;}; 
    virtual int T2_ctr_rel_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return 1;};
    virtual std::pair<int,int> ctr_rel_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return std::make_pair(-1,-1) ;};
    virtual std::pair<int,int> ctr_abs_pos(){ throw std::runtime_error("Not defined in CtrOp_base class!"); return std::make_pair(-1,-1) ;};;

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

    CtrOp_same_T(std::string T1name_in, std::string Tout_name_in, std::pair<int,int> ctr_abs_pos_in,
                 std::pair<int,int> ctr_rel_pos_in, std::string ctr_type_in ):CtrOp_base(Tout_name_in, ctr_type_in),
    T1name_(T1name_in), ctr_abs_pos_(ctr_abs_pos_in),  ctr_rel_pos_(ctr_rel_pos_in) {};

    ~CtrOp_same_T(){};

    std::string T1name() override { return T1name_ ;}
    std::pair<int,int> ctr_rel_pos() override { return ctr_abs_pos_;};
    std::pair<int,int> ctr_abs_pos() override { return ctr_rel_pos_;};
};      

template<class DType>
class TensorPart{
   public:
    std::string name;
    std::shared_ptr<std::vector<std::string>> idxs;
    std::shared_ptr<std::vector<std::string>> id_ranges;
    std::shared_ptr<DType> Tdata ;
  
    TensorPart(){};

    TensorPart(std::shared_ptr<std::vector<std::string>> idxs_in,
               std::shared_ptr<std::vector<std::string>> id_ranges_in,
               std::shared_ptr<DType> Tdata_in ) {
                  id_ranges = id_ranges_in; 
                  idxs = idxs_in;
                  Tdata = Tdata_in;
                  get_name();
                };

     virtual void get_name();

    ~TensorPart(){};
};

template<class DType>
class CtrTensorPart : public TensorPart<DType> /*, public: std::enable_shared_from_this<CtrTensorPart>*/ {
   public:
  
    using TP = TensorPart<DType> ;
    using TensorPart<DType>::idxs ;
    using TensorPart<DType>::id_ranges ;
    using TensorPart<DType>::name ;
   
    std::shared_ptr<std::vector<std::string>> full_idxs;
    std::shared_ptr<std::vector<std::string>> full_id_ranges;
    std::shared_ptr<std::vector<int>> unc_pos;
    std::shared_ptr<std::map<int,int>> unc_rel_pos;
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors; 
    std::shared_ptr<DType> CTdata ;
    std::shared_ptr<std::vector<std::string>> required_Tblocks;
    bool contracted; 
    bool survive_independently;
    

    CtrTensorPart(){};   
  
    CtrTensorPart(std::shared_ptr<std::vector<std::string>> full_idxs_in,
                  std::shared_ptr<std::vector<std::string>> full_id_ranges_in,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos_in,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors_in ) {

                  full_id_ranges = full_id_ranges_in; 
                  full_idxs = full_idxs_in;
                  ctrs_pos = ctrs_pos_in;
                  ReIm_factors = ReIm_factors_in;
                  contracted = false;
                  get_ctp_idxs_ranges();
                  get_name();
                  CTdata = std::make_shared<DType>(); 
                  required_Tblocks = std::make_shared<std::vector<std::string>>(0);
                };

     void get_ctp_idxs_ranges();

     void get_name() override;

     std::string myname(){ return name;};

     std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);

     void FullContract(std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap,
                       std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list_new );
    
    ~CtrTensorPart(){};
};


template<class DType>
class CtrMultiTensorPart : public TensorPart<DType> {
   public :

   using TensorPart<DType>::idxs ;
   using TensorPart<DType>::id_ranges ;
   using TensorPart<DType>::name ;

    std::string MTname;
    std::shared_ptr<std::vector<std::string>> full_idxs;
    std::shared_ptr<std::vector<std::string>> full_id_ranges;

    std::shared_ptr<std::vector<int>> Tsizes_cml;
    std::shared_ptr<std::vector<int>> all_unc_pos;
    std::shared_ptr<std::map<int,int>> all_unc_rel_pos;
    std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors;      
    std::shared_ptr<DType> CTdata ;
    std::shared_ptr<std::vector<std::string>> required_Tblocks;
    bool contracted;

    std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart<DType>>>> CTP_vec; 
    std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos;      

    CtrMultiTensorPart(){};

    CtrMultiTensorPart(std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart<DType>>>> CTP_vec_in, 
                       std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos_in ) {
                     
                       CTP_vec = CTP_vec_in;
                       cross_ctrs_pos = cross_ctrs_pos_in;
                       contracted = false;

                       full_idxs       = std::make_shared<std::vector<std::string>>(0);
                       full_id_ranges  = std::make_shared<std::vector<std::string>>(0);
                       Tsizes_cml      = std::make_shared<std::vector<int>>(0);
                       all_ctrs_pos    = std::make_shared<std::vector<std::pair<int,int>>>(0); 

                       int cml_size = 0;
                       for (auto ctp : *CTP_vec){
                         Tsizes_cml->push_back(cml_size);
                           
                         full_idxs->insert(full_idxs->end() , ctp->full_idxs->begin(), ctp->full_idxs->end());
                         
                         full_id_ranges->insert(full_id_ranges->end(),  ctp->full_id_ranges->begin(), ctp->full_id_ranges->end());

                         for (auto relctr : *ctp->ctrs_pos )
                           all_ctrs_pos->push_back( std::make_pair(relctr.first+Tsizes_cml->back(), relctr.second+Tsizes_cml->back()));
                         
                         cml_size+=ctp->full_idxs->size(); 
                         
                       }

                       for (auto cctr : *cross_ctrs_pos_in){
                         all_ctrs_pos->push_back(std::make_pair(Tsizes_cml->at(cctr.first.first)+cctr.first.second, Tsizes_cml->at(cctr.second.first)+cctr.second.second));
                       }
                       ReIm_factors = std::make_shared<std::vector<std::pair<int,int>>>(all_ctrs_pos->size());
                       get_cmtp_idxs_ranges();
                       get_name();
                     };
                                
    std::string myname(){ return name;};
    // NOTE : the order of the contractions in the name from get_name does _not_ match the ordering of the tensors. 
    // This greatly simplifies mapping from the gammas to the A-tensors
    void get_name() override;

    void get_cmtp_idxs_ranges();

    std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);

    void FullContract(std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap,
                      std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list_new );

    std::shared_ptr<CtrMultiTensorPart<DType>>
    Binary_Contract_diff_tensors_MT( std::string T1, std::string T2, std::pair<int,int> ctr_todo,
                                     std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>>> > Tmap,
                                     std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list_new );

    std::shared_ptr<CtrTensorPart<DType>>
    Binary_Contract_diff_tensors( std::pair<std::pair<int,int>, std::pair<int,int>> cross_ctr,
                                  std::pair<int,int>  ctr_todo,
                                  std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap,
                                  std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base> >> ACompute_list);

};

#endif
