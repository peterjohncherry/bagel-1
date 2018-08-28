#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename  DataType>
GammaInfo<DataType>::GammaInfo ( shared_ptr<CIVecInfo_Base> Bra_info, shared_ptr<CIVecInfo_Base> Ket_info,
                                 const vector<bool>& full_aops_vec, const vector<string>& full_idx_ranges,
                                 vector<int>& idxs_pos, shared_ptr<map<string, shared_ptr<GammaInfo_Base>>>& Gamma_map ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_GAMMA_INFO
 cout << "GammaInfo_Base:GammaInfo ; raw vector version" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  order_ = idxs_pos.size();
  Bra_info_ = Bra_info; 
  Ket_info_ = Ket_info; 

  id_ranges_ = make_shared<vector<string>>(order_);
  aops_      = make_shared<vector<bool>>(order_);

  {
  vector<string>::iterator ir_it = id_ranges_->begin();
  auto a_it = aops_->begin();
  for ( vector<int>::iterator ip_it = idxs_pos.begin(); ip_it != idxs_pos.end() ; ++ip_it, ++a_it, ++ir_it  ) {
    *ir_it = full_idx_ranges[*ip_it];
    *a_it  = full_aops_vec[*ip_it];
  }
  }

  sigma_id_ranges_ = make_shared<vector<string>>(order_+1);
  sigma_id_ranges_->front()= Bra_info_->name();
  copy( id_ranges_->begin(), id_ranges_->end(), (sigma_id_ranges_->begin()+1) );
  
  name_ = get_gamma_name( full_idx_ranges, full_aops_vec, idxs_pos, Bra_info_->name(), Ket_info_->name() ); 
  sigma_name_ = "S_"+name_;

  if (  (order_ > 2 ) && ( order_ % 2 == 0 ) ) {

    prev_gammas_ = vector<string>(order_/2);
    prev_sigmas_ = vector<string>(order_/2);

    for (int ii = 2 ; ii != order_; ii+=2 ){

      shared_ptr<vector<int>> prev_gamma_idxs_pos = make_shared<vector<int>> ( idxs_pos.begin()+ii, idxs_pos.end());
      prev_gammas_[(ii/2)-1] = get_gamma_name( full_idx_ranges, full_aops_vec, *prev_gamma_idxs_pos, Bra_info_->name(), Ket_info_->name());
      prev_sigmas_[(ii/2)-1] = "S_"+prev_gammas_[(ii/2)-1];

      if ( Gamma_map->find( prev_gammas_[(ii/2)-1] )  == Gamma_map->end() ){ //TODO fix Bra_info for rel case
        shared_ptr<GammaInfo<DataType>> prev_gamma_ = make_shared<GammaInfo<DataType>>( Bra_info_, Ket_info_, full_aops_vec, full_idx_ranges, *prev_gamma_idxs_pos, Gamma_map );
        Gamma_map->emplace( prev_gammas_[(ii/2)-1], prev_gamma_ );
      }
      if(  ii == 2 )
        prev_Bra_info_ = Gamma_map->at( prev_gammas_[0] )->Bra_info() ;
    }
  }

}
//////////////////////////////////////////////////////////////////////////
template class GammaInfo<double>;
template class GammaInfo<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
