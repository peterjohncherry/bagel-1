#include <bagel_config.h>
#include <map>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

using namespace std;
using namespace WickUtils;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO extend for rel case (Bra and Ket can vary,  prev gammas and
// prev_Bra_info should be constructed from idxs_pos, full_idxs_ranges, full_aops and Ket.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename  DataType>
GammaInfo<DataType>::GammaInfo ( shared_ptr<CIVecInfo<DataType>> Bra_info, shared_ptr<CIVecInfo<DataType>> Ket_info,
                                 shared_ptr<const vector<bool>> full_aops_vec, shared_ptr<const vector<string>> full_idx_ranges,
                                 shared_ptr<vector<int>> idxs_pos  ,
                                 shared_ptr<map<string, shared_ptr<GammaInfo<DataType>>>>& Gamma_map ) :
                                 order_(idxs_pos->size()), Bra_info_(Bra_info), Ket_info_(Ket_info) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// cout << "GammaInfo<DataType>:GammaInfo" << endl;
  id_ranges_ = make_shared<vector<string>>(idxs_pos->size());
  aops_      = make_shared<vector<bool>>(idxs_pos->size());
  for (int ii = 0 ; ii != idxs_pos->size(); ii++ ){
    id_ranges_->at(ii) = full_idx_ranges->at(idxs_pos->at(ii));
    aops_->at(ii)      = full_aops_vec->at(idxs_pos->at(ii));
  }

  sigma_id_ranges_ = make_shared<vector<string>>(id_ranges_->size()+1);
  sigma_id_ranges_->at(0) = Bra_info_->name();
  for ( int  ii = 1 ; ii != sigma_id_ranges_->size() ; ii++ )
    sigma_id_ranges_->at(ii) = id_ranges_->at(ii-1);

  name_ = get_gamma_name( full_idx_ranges, full_aops_vec, idxs_pos, Bra_info_->name(), Ket_info_->name() ); 
  sigma_name_ = "S_"+name_;

  if ( (idxs_pos->size() > 2 ) && ( idxs_pos->size() % 2 == 0 ) ) {

    prev_gammas_ = vector<string>(idxs_pos->size()/2);
    prev_sigmas_ = vector<string>(idxs_pos->size()/2);

    for (int ii = 2 ; ii != idxs_pos->size() ; ii+=2 ){

      shared_ptr<vector<int>> prev_gamma_idxs_pos = make_shared<vector<int>> ( idxs_pos->begin()+ii, idxs_pos->end());
      prev_gammas_[(ii/2)-1] = get_gamma_name( full_idx_ranges, full_aops_vec, prev_gamma_idxs_pos, Bra_info_->name(), Ket_info_->name());
      prev_sigmas_[(ii/2)-1] = "S_"+prev_gammas_[(ii/2)-1];

      if ( Gamma_map->find( prev_gammas_[(ii/2)-1] )  == Gamma_map->end() ){ //TODO fix Bra_info for rel case
        shared_ptr<GammaInfo<DataType>> prev_gamma_ = make_shared<GammaInfo<DataType>>( Bra_info_, Ket_info_, full_aops_vec, full_idx_ranges, prev_gamma_idxs_pos, Gamma_map );
        Gamma_map->emplace( prev_gammas_[(ii/2)-1], prev_gamma_ );
      }

      if(  ii == 2 )
        prev_Bra_info_ = Gamma_map->at( prev_gammas_[0] )->Bra_info() ;

    }
    print_vector( prev_gammas_ , "prev gammas of " + name_ ); cout << endl;
  }

}

//////////////////////////////////////////////////////////////////////////
template class GammaInfo<double>;
template class GammaInfo<std::complex<double>>;
/////////////////////////////////////////////////////////////////////////
