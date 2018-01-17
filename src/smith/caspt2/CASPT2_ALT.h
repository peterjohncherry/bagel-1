//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.h
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __SRC_SMITH_CASPT2_ALT_H
#define __SRC_SMITH_CASPT2_ALT_H

#include <src/smith/caspt2/CASPT2.h>
#include <src/ci/fci/fci.h>
#include <src/smith/multitensor.h>
#include <src/smith/wicktool/system_info.h>
#include <src/smith/wicktool/expression_computer.h>

namespace bagel {
namespace SMITH { 

namespace CASPT2 { class CASPT2; };
namespace TensOp_Computer { class TensOp_Computer<double>; };

namespace CASPT2_ALT{

class CASPT2_ALT {
  
  public:

 
   void set_range_info(std::shared_ptr<std::vector<int>> states_of_interest ) ;
   void set_target_info( std::shared_ptr<std::vector<int>> states_of_interest) ;
   void calculate_term( std::vector<int>& Target_states_list, std::vector<std::pair<std::vector<std::string>,double>>& BK_info_list, std::string term_type ) ;
   void solve();
   void Set_Tensor_Ops_Data( std::string op_name, std::string bra_name , std::string ket_name ) ;

    CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in);
    ~CASPT2_ALT() {};
    

};

}//end CASPT2_ALT namespace

}//end SMITH namespace 
}//end bagel namespace
#endif

