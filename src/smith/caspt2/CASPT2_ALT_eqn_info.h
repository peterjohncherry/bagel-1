#ifndef __SRC_SMITH_CASPT2_ALT_EQN_INFO_H
#define __SRC_SMITH_CASPT2_ALT_EQN_INFO_H

#include<stdlib.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include<utility>
#include<tuple>
#include<string>
#include<memory>
#include<map>
#include<numeric>

namespace bagel {
namespace SMITH {

namespace CASPT2_ALT_EQN_INFO {

   std::string flip(std::string idx);
   std::shared_ptr<std::vector<std::string>> ijkl_to_klij(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> ijkl_to_jilk(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> ijkl_to_lkji(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> ijkl_to_ijlk_block(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> ijkl_to_jikl_block(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> ijkl_to_jilk_block(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> bbbb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> bbaa_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> aabb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
   std::shared_ptr<std::vector<std::string>> identity(std::shared_ptr<std::vector<std::string>> invec) ;
   bool NotAllAct(std::shared_ptr<std::vector<std::string>> ranges);
   bool always_true(std::shared_ptr<std::vector<std::string>> ranges);
   std::vector<std::tuple<std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> set_2el_symmfuncs();

}//end of CASPT2_EQN_INFO namespace:
}
}
#endif
