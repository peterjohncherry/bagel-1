#ifndef __SRC_PROP_PROPTOOL_gamma_info_H
#define __SRC_PROP_PROPTOOL_gamma_info_H

#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/states_info.h>
#include <src/prop/proptool/algebraic_manipulator/range_block_info.h>

using namespace WickUtils;

template< typename DataType> 
class GammaInfo {

   private :
     int order_ ;

     //Follows convention : <I | i j | K > < K | .........| J >
     std::shared_ptr<CIVecInfo<DataType>> Bra_info_;      // < I |
     std::shared_ptr<CIVecInfo<DataType>> Ket_info_;      // | J >
     std::shared_ptr<CIVecInfo<DataType>> prev_Bra_info_; // < K |

     std::shared_ptr<std::vector<std::string>> id_ranges_ ;
     std::shared_ptr<std::vector<std::string>> sigma_id_ranges_ ;

     std::shared_ptr<std::vector<bool>> aops_ ;
     std::pair<DataType, DataType> factor_ ;

     std::string name_;
     std::string sigma_name_;

     std::vector<std::string> prev_gammas_;
     std::vector<std::string> prev_sigmas_;

   public :
     GammaInfo( std::shared_ptr<CIVecInfo<DataType>> Bra_info, std::shared_ptr<CIVecInfo<DataType>> Ket_info,
                std::shared_ptr<const std::vector<bool>> full_aops_vec, std::shared_ptr<const std::vector<std::string>> full_idx_ranges,
                std::shared_ptr<std::vector<int>> idxs_pos, std::shared_ptr<std::map< std::string, std::shared_ptr<GammaInfo<DataType>>>>& Gamma_map_in );
     GammaInfo(){};
    ~GammaInfo(){};

     int order() { return order_ ; };

     std::shared_ptr<std::vector<std::string>> id_ranges(){ return id_ranges_ ;};
     std::shared_ptr<std::vector<std::string>> sigma_id_ranges(){ return sigma_id_ranges_ ;};

     std::string name() { return name_; };
     std::string sigma_name() { return sigma_name_; };

     std::vector<std::string> prev_sigmas() { return prev_sigmas_; };
     std::string prev_sigmas(int ii ) { return prev_sigmas_[ii]; };
     std::string prev_sigma_name() { return prev_sigmas_[0]; };

     std::vector<std::string> prev_gammas() { return prev_gammas_; };
     std::string prev_gammas(int ii ) { return prev_gammas_[ii] ; };
     std::string prev_gamma_name() { return prev_gammas_[0]; };

     std::string Bra_name() { return Bra_info_->name(); };
     std::string Ket_name() { return Ket_info_->name(); };
     std::string Prev_Bra_name() { return prev_Bra_info_->name(); };

     std::shared_ptr<CIVecInfo<DataType>> Bra_info() {return  Bra_info_; };
     std::shared_ptr<CIVecInfo<DataType>> prev_Bra_info(){ return prev_Bra_info_; };
     std::shared_ptr<CIVecInfo<DataType>> Ket_info() { return Ket_info_; };

     int Bra_nalpha() const {return  Bra_info_->nalpha(); };
     int prev_Bra_nalpha() const { return prev_Bra_info_->nalpha(); };
     int Ket_nalpha() const { return Ket_info_->nalpha(); };

     int Bra_nbeta() const {return  Bra_info_->nbeta(); };
     int prev_Bra_nbeta() const { return prev_Bra_info_->nbeta(); };
     int Ket_nbeta() const { return Ket_info_->nbeta(); };

};

#endif
