#include <src/smith/wicktool/WickUtils.h>
#ifndef __SRC_SMITH_CtrTensOp_H
#define __SRC_SMITH_CtrTensOp_H

//#include "WickUtils.h"


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
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors; 
    std::shared_ptr<DType> CTdata ;
    bool contracted; 

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
                  get_name();;
                  CTdata = std::make_shared<DType>(); 
                };

     void get_ctp_idxs_ranges();
     void get_name() override;
     std::string myname(){ return name;};
     std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);
     void FullContract(std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap, std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );
     std::shared_ptr<DType> Binary_Contract_same_tensor(std::pair<int,int> ctr_todo, std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );

     std::shared_ptr<CtrTensorPart<DType>> Binary_Contract_same_tensor( std::string T1name , std::pair<int,int> ctr_todo,
                                                                        std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap ,
                                                                        std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );

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
    std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos;      
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors;      
    std::shared_ptr<DType> CTdata ;
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
                         full_idxs->insert(full_idxs->begin() , ctp->full_idxs->begin(), ctp->full_idxs->end());
                         full_id_ranges->insert(full_id_ranges->begin(),  ctp->full_id_ranges->begin(), ctp->full_id_ranges->end());

                         for (auto relctr : *ctp->ctrs_pos )
                           all_ctrs_pos->push_back( std::make_pair(relctr.first+Tsizes_cml->back(), relctr.second+Tsizes_cml->back()));

                         cml_size+=ctp->full_idxs->size(); 
                         
                       }

                       for (auto cctr : *cross_ctrs_pos_in){
                         all_ctrs_pos->push_back(std::make_pair(Tsizes_cml->at(cctr.first.first)+cctr.first.second, Tsizes_cml->at(cctr.second.first)+cctr.second.second));
                       }
                       ReIm_factors = std::make_shared<std::vector<std::pair<int,int>>>(all_ctrs_pos->size());

                       get_name();
                       //std::cout << "MTname = " << name << std::endl << std::endl;
                     };
                                
    std::string myname(){ return name;};
    void get_name() override;
    void get_name_readable();

    std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);

    void FullContract(std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>> >> Tmap, std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );

    std::shared_ptr<CtrMultiTensorPart<DType>> Binary_Contract_diff_tensors_MT(std::string T1, std::string T2, std::pair<int,int> ctr_todo,
                                                                               std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>>> > Tmap,
                                                                               std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );

    std::shared_ptr<CtrTensorPart<DType>> Binary_Contract_diff_tensors(std::string T1, std::string T2, std::pair<int,int> ctr_todo,
                                                                               std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart<DType>>> > Tmap,
                                                                               std::shared_ptr<std::vector< std::tuple<std::string,std::string,std::pair<int,int>> >> Acompute_list );

};


#endif
