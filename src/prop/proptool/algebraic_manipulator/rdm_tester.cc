#include <bagel_config.h>
#include <src/prop/proptool/proputils.h>
#include <src/prop/proptool/algebraic_manipulator/rdm_tester.h>

using namespace std;
 
///////////////////////////////////////////////////////////////////////////////////////////////////////
void RDM_Tester::make_gamma_info( int order, int bra_num, int ket_num  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RDM_TESTER
cout << "RDM_Tester::make_gamma_info" << endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////

  string bra_name = states_info_->civec_names( bra_num )->front();
  string ket_name = states_info_->civec_names( ket_num )->front();
  
  make_gamma_info( order, bra_name, ket_name ); 

  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void RDM_Tester::make_gamma_info( int order, string bra_name, string ket_name  ) {
///////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_RDM_TESTER
cout << "RDM_Tester::make_gamma_info" << endl;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////

   vector<bool> gamma_aops(order*2);
   for ( auto ga_it = gamma_aops.begin(); ga_it != gamma_aops.end() ; ++ga_it  ) {
     *ga_it = true;
      ++ga_it;
      *ga_it = false;
   }

   vector<string> gamma_rngs( order, "a" );
   vector<int> idxs_pos( order );
   iota( idxs_pos.begin(), idxs_pos.end(), 0 );
    
   std::shared_ptr<CIVecInfo_Base> bra_info = states_info_->civec_info(bra_name);
   std::shared_ptr<CIVecInfo_Base> ket_info = states_info_->civec_info(ket_name);

   string gamma_name = get_gamma_name( gamma_rngs, gamma_aops, idxs_pos, bra_name, ket_name );

   if ( gamma_map_->find(gamma_name) != gamma_map_->end() ){
     auto new_gamma = make_shared<GammaInfo<double>>( bra_info, ket_info, gamma_aops, gamma_rngs, idxs_pos, gamma_map_);  
     gamma_map_->emplace( gamma_name, new_gamma ); 
   } 

   return;
}
