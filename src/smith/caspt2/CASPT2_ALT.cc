//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/caspt2/CASPT2_ALT.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::TensOp_Computer; 
using namespace bagel::SMITH::Gamma_Computer; 
using namespace Tensor_Arithmetic_Utils;

////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in ) { 
////////////////////////////////////////////////////////////////////
  cout << " CASPT2_ALT::CASPT2_ALT::CASPT2_ALT" << endl;
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_range_info(shared_ptr<vector<int>> states_of_interest ) {
//////////////////////////////////////////////////////////////////////////////////////////////
  cout << "CASPT2_ALT::CASPT2_ALT::set_range_info " << endl;
  return;

}
//////////////////////////////////////////////////////////////////////////////////////////////
// defining this so we don't need to put determinants class into gamma generator etc.
// At the moment this is quite silly; the determinant space is the same for every space,
// however, this is not necessarily the case. This function will need to do something 
// different in the relativistic case.
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_target_info( shared_ptr<vector<int>> states_of_interest) {
//////////////////////////////////////////////////////////////////////////////////////////////
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::calculate_term( vector<int>& Target_states_list, vector<pair<vector<string>,double>>& BK_info_list, string term_type ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::solve() {
////////////////////////////////////////////////////////////////////
cout <<  " CASPT2_ALT::CASPT2_ALT::solve() " << endl;

 return;
} 
/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data( string op_name, string bra_name , string ket_name )  { 
/////////////////////////////////////////////////////////////////////////////////
cout << "CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data() " << endl;
  return;
}
#endif
