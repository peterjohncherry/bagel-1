#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/spin_manager.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace WickUtils;

using pint_vec = std::vector<std::pair<int,int>>;
using pstr_vec = std::vector<std::pair<std::string,std::string>>;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Spin_Manager::set_spin_info() {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  auto prvec = [&ss](shared_ptr<vector<bool>> invec){ cout << "[ " ; for (auto elem : *invec) { cout << ss(elem) << " " ;} cout << "]" ;};

  spin_max = (nact_orb < nact_el) ? nact_orb : nact_el;
  spin_max = (nact_el - nact_orb) > 0 ?  (nact_el - nact_orb)  : 0;
  
  spin_sectors = spin_sector_gen(nact_el , nact_orb);  
  spin_paths = make_shared<map< pair<int, pair<int,int>> , shared_ptr<vector<shared_ptr<vector<bool>>>> >>();
  for (int ii = 0; ii != (nidxs/2)+1 ; ii++){
    auto all_sc = all_spin_combinations(ii);
    for (pair<int,int> spin_sec : *spin_sectors){
      bool newsector = true;
      for (auto spin_comb : *all_sc) {
        if(sector_constrained_spin_combs( spin_sec.first, spin_sec.second, spin_comb)){
          if(!newsector){
            spin_paths->at(make_pair(ii, spin_sec))->push_back(spin_comb);
          } else {
            auto sscombs_vec = make_shared<vector<shared_ptr<vector<bool>>>>(1,spin_comb); 
            spin_paths->emplace(make_pair(ii, spin_sec), sscombs_vec);
            newsector = false;
          }
        }
      }
    }     
  }
  return; 
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Spin_Manager::bfunc(int xx, shared_ptr<vector<bool>> qvec, int nel, int length, shared_ptr<vector<shared_ptr<vector<bool>>>> all_sc) {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  auto xvec = make_shared<vector<bool>>(*qvec);
  if (nel == 0){
    all_sc->push_back(xvec);
    return ;
  } else {
    nel -=1;
    for (int ii = xx; ii != length-nel; ii++){ 
      xvec->at(ii)=true;
      bfunc(ii+1, xvec, nel, length, all_sc);
      xvec->at(ii)=false;
    }
    nel= nel-1;
  }
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Generates all true/false combinations for a vector which is _half_ the length of the total number of electrons.
//Thish is then used to build spin sector constrained transition lists. Should probably rename
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<vector<bool>>>>  Spin_Manager::all_spin_combinations(int length){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };

  auto xvec = make_shared<vector<bool>>(length, false);
  auto all_sc = make_shared<vector<shared_ptr<vector<bool>>>>(0) ;
  for (int nel = 0 ; nel != length+1 ; nel++){
     bfunc( 0, xvec, nel, length, all_sc); 
  }

  auto svecs = make_shared<vector<shared_ptr<vector<bool>>>>(0) ;
  for (auto plus_vec : *all_sc){
    for (auto kill_vec : *all_sc){
      auto buffvec = make_shared<vector<bool>>();
      buffvec->insert(buffvec->end(), plus_vec->begin(), plus_vec->end() );
      buffvec->insert(buffvec->end(), kill_vec->begin(), kill_vec->end() );
      svecs->push_back(buffvec);
    }
  }
  all_sc->swap(*svecs);
  return all_sc ;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Checks that the spin transition pathway never leaves the active spin_sectors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Spin_Manager::sector_constrained_spin_combs(int nalpha, int nbeta , shared_ptr<vector<bool>> gamma_comb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  for (int ii = gamma_comb->size()-1; ii != -1; ii-=2 ){
    if (gamma_comb->at(ii)){
      nalpha++;
    } else {
      nbeta++;
    }
    if( (nalpha < 0) || (nbeta < 0) ){
      return false; 
    }

    if (gamma_comb->at(ii-1)){
      nalpha--;
    }else {
      nbeta--;
    }
    if ((nbeta > nact_orb) || (nalpha > nact_orb)){
      return false; 
    }
  }
  
  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Checks that the spin transition pathway never leaves the active spin_sectors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Spin_Manager::sector_constrained_spin_combs(int nalpha, int nbeta , shared_ptr<vector<string>> gamma_comb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };
  for (int ii = gamma_comb->size()-1; ii != -1; ii-=2 ){
    if (gamma_comb->at(ii)=="A"){
      nalpha++;
    } else {
      nbeta++;
    }
    if( (nalpha < 0) || (nbeta < 0) ){
      return false; 
    }

    if (gamma_comb->at(ii-1)=="A"){
      nalpha--;
    }else {
      nbeta--;
    }
    if ((nbeta > nact_orb) || (nalpha > nact_orb)){
      return false; 
    }
  }
  return true;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//generates all possible spin sectors, should not be needed in Bagel
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<int,int>>> Spin_Manager::spin_sector_gen(int nact_el, int nact_orb){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  int max_spin = nact_el > nact_orb ? nact_orb : nact_el;
  int min_spin = nact_el - nact_orb > 0 ? nact_el- nact_orb : 0;

  auto spin_sectors = make_shared<vector<pair<int,int>>>(0);
  for (int ii = 0 ; ii != max_spin-min_spin+1; ii++)
    spin_sectors->push_back(make_pair(min_spin+ii, max_spin-ii));

  return spin_sectors;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Temporary fix, should either decide on having bools or strings for spins
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Spin_Manager::get_string_spin_paths(){
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto ss = [](bool spin ) { return spin ?  "A" : "B" ; };

  spin_paths_str =  make_shared<map< pair< int, pair<int,int> > , shared_ptr<vector<shared_ptr<vector<string>>>> >>();

  for (auto mapit= spin_paths->begin(); mapit!=spin_paths->end(); mapit++){
   auto boolspinvecs = mapit->second;
   auto strspinvecs = make_shared<vector<shared_ptr<vector<string>>>>(0);
   for (auto vec : *boolspinvecs){ 
     auto strspins = make_shared<vector<string>>(0);
     for (auto bspin : *vec) {
       strspins->push_back(ss(bspin)); 
     }
     strspinvecs->push_back(strspins);
    }
   spin_paths_str->emplace(mapit->first, strspinvecs);
  }
  return;
}
