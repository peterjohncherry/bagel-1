#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/braket.h>

using namespace std;
using namespace WickUtils;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BraKet_Base::BraKet_Base( std::shared_ptr<Op_Info> multiop_info, 
                          std::pair<double, double> factor, int bra_num, int ket_num,  std::string type) :
                          multiop_info_(multiop_info), factor_(factor), bra_num_(bra_num), ket_num_(ket_num),
                          type_(type) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_BRAKET_BASE
cout << "BraKet_Base::BraKet_Base" << endl;
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if (type_[0] == 'c' ) {
    name_ = "c_{I}"; 
  } else {
    name_ = ""; 
  }
  name_ += "<" + to_string(bra_num)+ "| " + multiop_info->op_full_name_  + " |" + to_string(ket_num) + ">";
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BraKet_Base::print_gamma_Atensor_contractions(shared_ptr<map<string, shared_ptr< map<string, shared_ptr<AContribInfo_Base> >>>> G_to_A_map,
		                                        bool has_orb_exc ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_BRAKET_BASE
  cout <<  "BraKet_Base::print_gamma_Atensor_contractions()" << endl; 
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  for( auto map_it = G_to_A_map->begin() ; map_it != G_to_A_map->end(); map_it++){
  
    cout << "====================================================" << endl;
    cout << map_it->first << endl;
    cout << "====================================================" << endl;
  
    if ( map_it->second->size() == 0 ) {
  
      cout << endl << "      - No Contributions! -" <<  endl << endl;
  
    } else {
  
      for( auto A_map_it = map_it->second->begin() ; A_map_it != map_it->second->end();  A_map_it++){
        cout <<  A_map_it->first << "  "; cout.flush();
        for ( int qq = 0; qq !=A_map_it->second->id_orders().size(); qq++) {
          cout << "{ " ; print_vector( A_map_it->second->id_order(qq), "" ); 
          cout << " (" << A_map_it->second->factor(qq).first <<  "," <<  A_map_it->second->factor(qq).first<<  ") }" ; cout <<endl;
        }
        cout << endl;
      }
    }
  }
}
