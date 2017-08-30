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

#include <tuple>
#include <src/smith/wicktool/equation.h>
#include <src/smith/caspt2/CASPT2.h>

namespace bagel {
namespace SMITH { 

namespace CASPT2 { class CASPT2; };

namespace CASPT2_ALT{

class CASPT2_ALT {
  
  public:
    
    std::shared_ptr<CASPT2::CASPT2> orig_cpt2;
    std::shared_ptr<const SMITH_Info<double>> info_;

    CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in);
    std::shared_ptr<Dvec> cc_; 
    int  nelea_ ;
    int  neleb_ ;
    int  ncore_ ;
    int  norb_  ;
    int  nstate_;
    std::shared_ptr<Determinants> det_ ; 


    CASPT2_ALT(std::shared_ptr<const SMITH_Info<double>> ref_alt);

    ~CASPT2_ALT() {};
    
    void test();
    void compute_gamma12(const int MM, const int NN ) ;

    //std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
    void
    compute_gamma12_from_civec(std::shared_ptr<const Civec> cbra, std::shared_ptr<const Civec> cket) const ;

//  std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
    void
    compute_gamma12_last_step(std::shared_ptr<const Dvec> dbra, std::shared_ptr<const Dvec> dket, std::shared_ptr<const Civec> cibra) const ;

};
}//end CASPT2_ALT namespace


namespace WICKTOOLS{

class WICKTOOLS {
  
  public:
    
    std::shared_ptr<Dvec> cc_; 
    int  nelea_ ;
    int  neleb_ ;
    int  ncore_ ;
    int  norb_  ;
    int  nstate_;
    std::shared_ptr<Determinants> det_ ; 
    WICKTOOLS(std::shared_ptr<const SMITH_Info<double>> ref_alt);
    ~WICKTOOLS() {};
    
    void compute_gamma12(const int MM, const int NN ) ;

    //std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
    void
    compute_gamma12_from_civec(std::shared_ptr<const Civec> cbra, std::shared_ptr<const Civec> cket) const ;

//  std::tuple<std::shared_ptr<RDM<1>>, std::shared_ptr<RDM<2>>>
    void
    compute_gamma12_last_step(std::shared_ptr<const Dvec> dbra, std::shared_ptr<const Dvec> dket, std::shared_ptr<const Civec> cibra) const ;

};
}//end WICKTOOLS namespace


}//end SMITH namespace 
}//end bagel namespace
#endif

