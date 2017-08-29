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

#include <iostream>
#include <iomanip>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in ) { 
////////////////////////////////////////////////////////////////////
   
//   CASPT2::CASPT2 orig_cpt2_raw(orig_cpt2_in);
   orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
}
////////////////////////////////////////////////////////////////////
CASPT2_ALT::CASPT2_ALT::CASPT2_ALT(const double factor_in ) { 
////////////////////////////////////////////////////////////////////
    
    factor = factor_in;
}

/////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::test() { 
/////////////////////////////////
  auto weqn = make_shared<Equation<vector<double>>>();
  weqn->equation_build();
}

#endif
