#ifndef __SRC_PROP_PROPTOOL_CtrTensOp_H
#define __SRC_PROP_PROPTOOL_CtrTensOp_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/task_translator/tensor_algop_info.h>
#include <unordered_set>

class CtrTensorPart_Base  {
  public : //TODO change to private

    std::shared_ptr<std::vector<std::string>> full_idxs_;
    std::shared_ptr<std::vector<std::string>> full_id_ranges_;
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos_;
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors_;

    std::string name_;

    std::shared_ptr<std::vector<int>> unc_pos_;
    std::shared_ptr<std::map<int,int>> unc_rel_pos_;

    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_todo_;
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_done_;

    std::shared_ptr<std::vector<std::string>> unc_idxs_;
    std::shared_ptr<std::vector<std::string>> unc_id_ranges_;

    std::unordered_set<std::string> dependents_;
    std::unordered_set<std::string> dependencies_;
    bool got_data_;
    bool got_compute_list_;
    bool contracted_;
 
    int size_;

    std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>> CTP_vec_;
  public :
    CtrTensorPart_Base()  { //TODO Fix this rubbish 
                    full_idxs_      = std::make_shared< std::vector<std::string>>(0);
                    full_id_ranges_ = std::make_shared< std::vector<std::string>>(0);
                    ctrs_pos_       = std::make_shared< std::vector<std::pair<int,int>>>(0);
                    ctrs_todo_      = std::make_shared< std::vector<std::pair<int,int>>>(0);
                    ReIm_factors_   = std::make_shared< std::vector<std::pair<int,int>>>(ctrs_pos_->size());
                    contracted_     = false;
                    got_compute_list_ = false;
                   };


    CtrTensorPart_Base( std::shared_ptr<std::vector<std::string>> full_idxs,
                        std::shared_ptr<std::vector<std::string>> full_id_ranges,
                        std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos,
                        std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors ) :
                        full_id_ranges_(full_id_ranges), full_idxs_(full_idxs), ctrs_pos_(ctrs_pos),
                        ReIm_factors_(ReIm_factors), got_data_(false), size_( full_idxs_->size()),
                        name_(WickUtils::get_ctp_name(*full_idxs_, *full_id_ranges_, *ctrs_pos_ )) {  
                          ctrs_todo_ = std::make_shared<std::vector<std::pair<int,int>>>(*ctrs_pos);
                          ctrs_done_ = std::make_shared<std::vector<std::pair<int,int>>>(0);
                          got_data_ = false;
                          get_ctp_idxs_ranges();
                          CTP_vec_ = std::make_shared<std::vector<std::shared_ptr<CtrTensorPart_Base>>>(1);
                        } 

    ~CtrTensorPart_Base(){};

    std::string name() { return name_; }
    std::shared_ptr<std::vector<std::string>> full_idxs() { return  full_idxs_; }
    std::shared_ptr<std::vector<std::string>> full_id_ranges() { return full_id_ranges_; }
    std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors() { return ReIm_factors_; }

    std::shared_ptr<std::vector<int>> unc_pos() { return unc_pos_; }
    std::shared_ptr<std::vector<std::string>> unc_idxs() { return unc_idxs_ ; }
    std::shared_ptr<std::vector<std::string>> unc_id_ranges() { return unc_id_ranges_; }
    std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos() { return ctrs_pos_; }

    std::string full_idxs(int ii) { return  full_idxs_->at(ii); }
    std::string full_id_ranges(int ii) { return full_id_ranges_->at(ii); }

    std::shared_ptr<const std::vector<int>> unc_pos() const { return unc_pos_; }
    std::shared_ptr<std::vector<std::string>> unc_id_ranges() const { return unc_id_ranges_; }

    int unc_pos( int ii ) const { return unc_pos_->at(ii); }
    std::string unc_id_ranges( int ii ) const { return unc_id_ranges_->at(ii) ; }

    bool got_compute_list(){ return got_compute_list_; }
    void got_compute_list( bool val ){ got_compute_list_ = val; }
    
    std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>> CTP_vec() { return CTP_vec_; }
    
    int size() { return full_idxs_->size() ; } 
  
    void get_name();
    void get_ctp_idxs_ranges();

    std::string get_next_name(std::shared_ptr<std::vector<std::pair<int,int>>> new_ctrs_pos);

    virtual void FullContract( std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart_Base>>> CTP_map,
                               std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list ,
                               std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map) = 0 ; 
//                              { throw std::logic_error(" Calling full contract from CtrTensorPart_Base !! Aborting !! ");  } ; 
    
    virtual std::pair<int,int> get_pre_contract_ctr_rel_pos( std::pair<int,int>& ctr_pos ) = 0;
};

template<typename DataType>
class CtrTensorPart : public  CtrTensorPart_Base , public std::enable_shared_from_this<CtrTensorPart<DataType>>   {
   public:

    CtrTensorPart() : CtrTensorPart_Base() {} 

    CtrTensorPart(std::shared_ptr<std::vector<std::string>> full_idxs,
                  std::shared_ptr<std::vector<std::string>> full_id_ranges,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ctrs_pos,
                  std::shared_ptr<std::vector<std::pair<int,int>>> ReIm_factors ) :
                  CtrTensorPart_Base( full_idxs, full_id_ranges, ctrs_pos, ReIm_factors) {};

     void FullContract( std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart_Base> >> CTP_map,
                        std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list ,
                        std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> ACompute_map );

    std::pair<int,int> get_pre_contract_ctr_rel_pos( std::pair<int,int>& ctr_pos );

    ~CtrTensorPart(){};
};


template<typename DataType>
class CtrMultiTensorPart : public CtrTensorPart_Base , public  std::enable_shared_from_this<CtrMultiTensorPart<DataType>>  {
   public :

    std::shared_ptr<std::vector<int>> Tsizes_cml;
    std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos_;

    bool get_compute_list_from_reordered_tens_;
    std::string  reordered_tens_name_;
    std::shared_ptr<std::vector<int>> reordering_; 

    CtrMultiTensorPart(){};

    CtrMultiTensorPart( std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>> CTP_vec,
                        std::shared_ptr<std::vector<std::pair<std::pair<int,int>, std::pair<int,int>> >> cross_ctrs_pos_in  )
                        : CtrTensorPart_Base() {
                         std::cout << "CtrMultiTensorPart::CtrMultiTensorPart" << std::endl; 
                         int counter = 0;
                         CTP_vec_         = CTP_vec;
                         cross_ctrs_pos_ = cross_ctrs_pos_in;
                         Tsizes_cml      = std::make_shared<std::vector<int>>(0);
                         ctrs_pos_       = std::make_shared<std::vector<std::pair<int,int>>>(0);
                         got_compute_list_ = false;
                         get_compute_list_from_reordered_tens_ = false;

                         int cml_size = 0;
                         for (std::shared_ptr<CtrTensorPart_Base> ctp : *CTP_vec_){
                           Tsizes_cml->push_back(cml_size);
                           full_idxs_->insert(full_idxs_->end() , ctp->full_idxs_->begin(), ctp->full_idxs_->end());
                           full_id_ranges_->insert(full_id_ranges_->end(),  ctp->full_id_ranges_->begin(), ctp->full_id_ranges_->end());
                           std::cout << "ctp->name() = " <<  ctp->name() << std::endl;
                           WickUtils::print_pair_vector( *ctp->ctrs_pos_,  ctp->name() + " ctrs_pos in constructor"  ); std::cout << std::endl; 
                           WickUtils::print_vector( *ctp->full_idxs_, ctp->name() + " full_idxs in constructor"  ); std::cout << std::endl; 
                           for (auto relctr : *ctp->ctrs_pos_ )
                             ctrs_pos_->push_back( std::make_pair(relctr.first+Tsizes_cml->back(), relctr.second+Tsizes_cml->back()));

                           WickUtils::print_pair_vector( *ctp->ctrs_pos_, ctp->name() +  " ctrs_pos in constructor"  ); std::cout << std::endl; 
                           cml_size+=ctp->full_idxs_->size();
                           std::cout << "cml_size = " << cml_size << std::endl;
                         }
                       
          
                         WickUtils::print_vector( *Tsizes_cml, "Tsizes_cml" ); std::cout << std::endl; 
                         WickUtils::print_pair_pair_vector( *cross_ctrs_pos_in, "cross_ctrs_pos_in" ); std::cout << std::endl;
                         for (auto cctr : *cross_ctrs_pos_in) { 
                           ctrs_pos_->push_back(std::make_pair(Tsizes_cml->at(cctr.first.first)+cctr.first.second, Tsizes_cml->at(cctr.second.first)+cctr.second.second));
                         }
                         WickUtils::print_pair_vector( *ctrs_pos_, " full ctrs_pos in constructor"  ); std::cout << std::endl; 

                         get_ctp_idxs_ranges();
                         this->name_ = WickUtils::get_ctp_name(*full_idxs_, *full_id_ranges_, *ctrs_pos_ ) ;
                         std::cout << "finished building cmtp" << std::endl;
                       };


    //fix names; CTP_vec should at least be private
    std::shared_ptr<std::vector<std::shared_ptr<CtrTensorPart_Base>>> get_CTP_vec() const { return CTP_vec_; }  ;

    std::pair<int,int> get_pre_contract_ctr_rel_pos( std::pair<int,int>& ctr_pos ) { throw std::logic_error("not implemented yet!! Aborting!!" );  return std::make_pair(-1, -1 ); };

    void use_new_order_compute_list( std::shared_ptr<std::vector<int>> reordering, std::string reordered_tens_name ) {
             get_compute_list_from_reordered_tens_ = true;
             reordering_ = reordering;
             reordered_tens_name_ = reordered_tens_name; return; } 
                                                                                         

    void FullContract( std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart_Base> >> Tmap,
                       std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list ,
                       std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);

    std::shared_ptr<CtrTensorPart<DataType>>
    Binary_Contract_diff_tensors( std::pair<std::pair<int,int>, std::pair<int,int>> cross_ctr,
                                  std::pair<int,int>  ctr_todo,
                                  std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart_Base>>> CTP_map,
                                  std::shared_ptr<std::vector<std::shared_ptr<CtrOp_base> >> ACompute_list,
                                  std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);

    // TODO : Remove, this should not be needed once the mutiltensor nesting is done correctly.
    std::shared_ptr<CtrMultiTensorPart<DataType>>
    Binary_Contract_diff_tensors_MT( std::string T1, std::string T2, std::pair<int,int> ctr_todo,
                                     std::shared_ptr<std::map<std::string,std::shared_ptr<CtrTensorPart_Base>>> CTP_map,
                                     std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >> Acompute_list_new ,
                                     std::shared_ptr<std::map<std::string, std::shared_ptr<std::vector< std::shared_ptr<CtrOp_base> >>>> Acompute_map);


};

#endif
