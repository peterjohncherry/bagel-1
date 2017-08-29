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

#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/wicktool/equation.h>
#include <src/smith/caspt2/CASPT2.h>

//namespace bagel {
//namespace SMITH { 

//namespace CASPT2_ALT{

class CASPT2_ALT {
  
  public:
//  const CASPT2::CASPT2& orig_cas;

    //CASPT2_ALT(){};
    CASPT2_ALT(/*const CASPT2::CASPT2& orig_cas_in*/ double xx);

    ~CASPT2_ALT() {};

    void test();

};
//}
//}
//}

//}
#endif

