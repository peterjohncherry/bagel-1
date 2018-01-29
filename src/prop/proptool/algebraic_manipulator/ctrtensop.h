#ifndef __SRC_PROP_PROPTOOL_CtrTensOp_H
#define __SRC_PROP_PROPTOOL_CtrTensOp_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/task_translation/tensor_algop_info.h>
#include <unordered_set>

template<typename DataType>
class CtrTensorPart_Base  {
  private :
    const std::string name_;
    const std::shared_ptr<const std::vector<std::string>> full_idxs_;                                                                              
    const std::shared_ptr<const std::vector<std::string>> full_id_ranges_;
    const std::pair<double,double> ReIm_factors_; 

    const std::vector<int> unc_pos_;
    const std::map<int,int> unc_rel_pos_;
    const std::vector<std::string> unc_idxs_;
    const std::vector<std::string> unc_id_ranges_;
    const std::vector<std::pair<int,int>> ctrs_pos_;      

    const std::shared_ptr<const std::vector<int>> unc_pos_ptr_;
    const std::shared_ptr<const std::vector<std::string>> unc_idxs_ptr_;
    const std::shared_ptr<const std::vector<std::string>> unc_id_ranges_ptr_;
    const std::shared_ptr<const std::vector<std::pair<int,int>>> ctrs_pos_ptr_;      

  public :
    CtrTensorPart_Base(std::string name,
                       std::shared_ptr<const std::vector<std::string>> full_idxs,                                                                              
                       std::shared_ptr<const std::vector<std::string>> full_id_ranges,
                       const std::pair<double,double>& ReIm_factors, 
                       std::vector<int>& unc_pos,
                       std::map<int,int>& unc_rel_pos,
                       std::vector<std::string>& unc_idxs,
                       std::vector<std::string>& unc_id_ranges,
                       std::vector<std::pair<int,int>>& ctrs_pos) :
                       name_(name), full_idxs_(full_idxs), full_id_ranges_(full_id_ranges), ReIm_factors_(ReIm_factors),
                       unc_pos_(unc_pos), unc_rel_pos_(unc_rel_pos), unc_idxs_(unc_idxs), unc_id_ranges_(unc_id_ranges),
                       ctrs_pos_(ctrs_pos),      
                       unc_pos_ptr_(std::make_shared<const std::vector<int>>(unc_pos)),
                       unc_idxs_ptr_(std::make_shared<const std::vector<std::string>>(unc_idxs)),
                       unc_id_ranges_ptr_(std::make_shared<const std::vector<std::string>>(unc_id_ranges)),
                       ctrs_pos_ptr_(std::make_shared<const std::vector<std::pair<int,int>>>(ctrs_pos)){}      
    ~CtrTensorPart_Base(){}; 

    std::shared_ptr<const std::vector<int>> unc_pos() const { return unc_pos_ptr_; };
    std::shared_ptr<const std::vector<std::string>> unc_idxs() const { return unc_idxs_ptr_ ; };
    std::shared_ptr<const std::vector<std::string>> unc_id_ranges() const { return unc_id_ranges_ptr_; };
    std::shared_ptr<const std::vector<std::pair<int,int>>> ctrs_pos() const { return ctrs_pos_ptr_; };      

};

template<typename DataType>
class CtrTensorPart /*, public: std::enable_shared_from_this<CtrTensorPart>*/ {
   public:
  
    std::string name;
    std::shared_ptr<std::vector<std::string>> unc_idxs;
    std::shared_ptr<std::vector<std::string>> unc_id_ranges;

    std::shared_ptr<std::vector<std::string>> full_idxs;
    std::shared_ptr<std::vector<std::string>> full_id_ranges;
    std::shared_ptr<std::vector<int>> unc_pos;
    std::shared_ptr<std::map<int,int>> unc_rel_pos;
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos;      

    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_todo;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_done;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors; 

    std::unordered_set<std::string> dependents;
    std::unordered_set<std::string> dependencies;
    bool got_data; 
    bool got_compute_list; 
    bool contracted; 
    
    CtrTensorPart(){ //TODO Fix this rubbish so it CMTP  uses the base class
                    full_idxs      = std::make_shared< std::vector<std::string>>(0);
                    full_id_ranges = std::make_shared< std::vector<std::string>>(0);
                    ctrs_pos       = std::make_shared< std::vector<std::pair<int,int>>>(0); 
                    ctrs_todo      = std::make_shared< std::vector<std::pair<int,int>>>(0);
                    ReIm_factors   = std::make_shared< std::vector<std::pair<int,int>>>(ctrs_pos->size());
                    contracted     = false;
                    got_compute_list = false;
                   };   
  
    CtrTensorPart(std::shared_ptr<std::vector<std::string>> full_idxs_in,
                  std::shared_ptr<std::vector<std::string>> full_id_ranges_in,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos_in,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors_in ) :
                  full_id_ranges(full_id_ranges_in), full_idxs(full_idxs_in), ctrs_pos(ctrs_pos_in),
                  ReIm_factors(ReIm_factors_in) {
                  ctrs_todo = std::make_shared<std::vector<std::pair<int,int>>>(*ctrs_pos_in);
                  ctrs_done = std::make_shared<std::vector<std::pair<int,int>>>(0);
                  got_data = false;
                  get_ctp_idxs_ranges();
                  get_name();
                };


     void get_ctp_idxs_ranges();

     void get_name();

     std::string myname(){ return name;};

     std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);

     void FullContract(std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DataType>> >> Tmap,
                       std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);
    
    std::pair<int,int> get_pre_contract_ctr_rel_pos( std::pair<int,int>& ctr_pos );

    ~CtrTensorPart(){};
};


template<typename DataType>
class CtrMultiTensorPart : public CtrTensorPart<DataType>  {
   public :

    //should be a better way....
    using CtrTensorPart<DataType>::name;
    using CtrTensorPart<DataType>::unc_idxs;
    using CtrTensorPart<DataType>::unc_id_ranges;
    using CtrTensorPart<DataType>::full_idxs;
    using CtrTensorPart<DataType>::full_id_ranges;
    using CtrTensorPart<DataType>::unc_pos;
    using CtrTensorPart<DataType>::unc_rel_pos;
    using CtrTensorPart<DataType>::ctrs_pos;      
    using CtrTensorPart<DataType>::ctrs_todo; 
    using CtrTensorPart<DataType>::ctrs_done; 
    using CtrTensorPart<DataType>::ReIm_factors; 
    using CtrTensorPart<DataType>::dependents;
    using CtrTensorPart<DataType>::dependencies;
    using CtrTensorPart<DataType>::got_data; 
    using CtrTensorPart<DataType>::got_compute_list; 
    using CtrTensorPart<DataType>::contracted;

    std::shared_ptr<std::vector<int>> Tsizes_cml;
    std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart<DataType>>>> CTP_vec; 
    std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos;      

    CtrMultiTensorPart(){};

    CtrMultiTensorPart( std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart<DataType>>>> CTP_vec_in, 
                        std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos_in  ) 
                        : CtrTensorPart<DataType>() {
                    
                         CTP_vec        = CTP_vec_in;
                         cross_ctrs_pos = cross_ctrs_pos_in;
                         Tsizes_cml     = std::make_shared<std::vector<int>>(0);
                         ctrs_pos       = std::make_shared<std::vector<std::pair<int,int>>>(0); 
                         got_compute_list = false;                         
 
                         int cml_size = 0;
                         for (std::shared_ptr<CtrTensorPart<DataType>> ctp : *CTP_vec){
                           Tsizes_cml->push_back(cml_size);
                           full_idxs->insert(full_idxs->end() , ctp->full_idxs->begin(), ctp->full_idxs->end());
                           full_id_ranges->insert(full_id_ranges->end(),  ctp->full_id_ranges->begin(), ctp->full_id_ranges->end());
                           for (auto relctr : *ctp->ctrs_pos )
                             ctrs_pos->push_back( std::make_pair(relctr.first+Tsizes_cml->back(), relctr.second+Tsizes_cml->back()));
                           cml_size+=ctp->full_idxs->size(); 
                         }
                         
                         for (auto cctr : *cross_ctrs_pos_in){
                           ctrs_pos->push_back(std::make_pair(Tsizes_cml->at(cctr.first.first)+cctr.first.second, Tsizes_cml->at(cctr.second.first)+cctr.second.second));
                         }
                         get_ctp_idxs_ranges();
                         get_name();
                       };
 
    void get_name(){ CtrTensorPart<DataType>::get_name(); return;};

    std::string myname(){ return CtrTensorPart<DataType>::myname(); };

    std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos) {return CtrTensorPart<DataType>::get_next_name(new_ctrs_pos); };

    void get_ctp_idxs_ranges() { CtrTensorPart<DataType>::get_ctp_idxs_ranges(); return;};

    void FullContract( std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DataType>> >> Tmap,
                       std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);

    std::shared_ptr<CtrTensorPart<DataType>>
    Binary_Contract_diff_tensors( std::pair<std::pair<int,int>, std::pair<int,int>> cross_ctr,
                                  std::pair<int,int>  ctr_todo,
                                  std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DataType>> >> Tmap,
                                  std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base> >> ACompute_list,
                                  std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);

    std::shared_ptr<CtrMultiTensorPart<DataType>>
    Binary_Contract_diff_tensors_MT( std::string T1, std::string T2, std::pair<int,int> ctr_todo,
                                     std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DataType>>> > Tmap,
                                     std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list_new ,
                                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);


};

#endif
