#include <bagel_config.h>
#ifdef COMPILE_SMITH
//#include <src/smith/caspt2/CASPT2_ALT_eqn_info.h>
#include <src/smith/caspt2/CASPT2_ALT.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

/////////////////////////////////////////////////////////////
std::string CASPT2_ALT::CASPT2_ALT::flip(std::string idx){
/////////////////////////////////////////////////////////////
  if(idx.back() == 'A') idx.back()='A';
  if(idx.back() == 'B') idx.back()='B';
return idx;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//BEWARE!!!: all symm ops use physicists' ordering of integrals; most publications define the ordering of the
//           T and lambda indexes using chemists' ordering.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_klij(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {invec->at(2), invec->at(3), invec->at(0), invec->at(1)} ;
  auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
  return transformed_vec_ptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_jilk(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {invec->at(1), invec->at(0), invec->at(3), invec->at(2)} ;
  auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
  return transformed_vec_ptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_lkji(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {invec->at(3), invec->at(2), invec->at(1), invec->at(0)} ;
  auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
  return transformed_vec_ptr;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_ijlk_block(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if ( invec->at(2) != invec->at(3) ) {
     std::vector<std::string> transformed_vec = {invec->at(0), invec->at(1), invec->at(3), invec->at(2)} ;
     auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
     return transformed_vec_ptr;
  } 
  return invec;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_jikl_block(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (  invec->at(0) != invec->at(1) ) {
     std::vector<std::string> transformed_vec = {invec->at(1), invec->at(0), invec->at(2), invec->at(3)} ;
     auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
     return transformed_vec_ptr;
  } 
  return invec;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::ijkl_to_jilk_block(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if ( (invec->at(2) != invec->at(3)) && ( invec->at(0) != invec->at(1) ) ){
    std::vector<std::string> transformed_vec = {invec->at(1), invec->at(0), invec->at(3), invec->at(2)} ;
    auto transformed_vec_ptr = std::make_shared<std::vector<std::string>>(transformed_vec);
    return transformed_vec_ptr;
  } 
  return invec;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<string>> CASPT2_ALT::CASPT2_ALT::bbbb_to_aaaa(shared_ptr<vector<string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {flip(invec->at(0)), flip(invec->at(1)), flip(invec->at(2)), flip(invec->at(3))} ;
  auto transformed_vec_ptr = make_shared<vector<string>>(transformed_vec);
  return transformed_vec_ptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<string>> CASPT2_ALT::CASPT2_ALT::bbaa_to_aaaa(shared_ptr<vector<string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {flip(invec->at(0)), flip(invec->at(1)), invec->at(2), invec->at(3)} ;
  auto transformed_vec_ptr = make_shared<vector<string>>(transformed_vec);
  return transformed_vec_ptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<string>> CASPT2_ALT::CASPT2_ALT::aabb_to_aaaa(shared_ptr<vector<string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<std::string> transformed_vec = {invec->at(0), invec->at(1), flip(invec->at(2)), flip(invec->at(3))} ;
  auto transformed_vec_ptr = make_shared<vector<string>>(transformed_vec);
  return transformed_vec_ptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<std::vector<std::string>> CASPT2_ALT::CASPT2_ALT::identity(std::shared_ptr<std::vector<std::string>> invec) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return invec;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CASPT2_ALT::CASPT2_ALT::NotAllAct(shared_ptr<vector<string>> ranges){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (auto  elem : *ranges) 
    if ( elem != "act")
      return true;
  return false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool CASPT2_ALT::CASPT2_ALT::always_true(shared_ptr<vector<string>> ranges){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >>   CASPT2_ALT::CASPT2_ALT::set_2el_symmfuncs(){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector< tuple< shared_ptr<vector<string>>(*)(shared_ptr<vector<string>>),int,int >> Two_el_symmfuncs;
  int one = 1;
  int mone = -1;
//  shared_ptr<vector<string>> (*symmop1)(shared_ptr<vector<string>>)  = &ijkl_to_klij;
//  shared_ptr<vector<string>> (*symmop2)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::ijkl_to_jilk;
//  shared_ptr<vector<string>> (*symmop3)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::ijkl_to_lkji;
//  shared_ptr<vector<string>> (*symmop4)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::ijkl_to_ijlk_block;
//  shared_ptr<vector<string>> (*symmop5)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::ijkl_to_jikl_block;
//  shared_ptr<vector<string>> (*symmop6)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::ijkl_to_jilk_block;
//  shared_ptr<vector<string>> (*symmop7)(shared_ptr<vector<string>>)  = &bagel::SMITH::CASPT2_ALT::CASPT2_ALT::identity;
//  Two_el_symmfuncs.push_back(tie(symmop1, one, mone));
//  Two_el_symmfuncs.push_back(tie(symmop2, one, one));
//  Two_el_symmfuncs.push_back(tie(symmop3, one, mone));
//  Two_el_symmfuncs.push_back(tie(symmop4, mone, mone));
//  Two_el_symmfuncs.push_back(tie(symmop5, mone, mone));
//  Two_el_symmfuncs.push_back(tie(symmop6, mone, mone));
//  Two_el_symmfuncs.push_back(tie(symmop7, one, one));

  return Two_el_symmfuncs;  
} 

#endif
