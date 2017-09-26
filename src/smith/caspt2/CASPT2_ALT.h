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
#include <src/ci/fci/fci.h>
#include <src/smith/wicktool/equation_tools.h>
#include <src/smith/multitensor.h>

namespace bagel {
namespace SMITH { 

namespace CASPT2 { class CASPT2; };
namespace Equation_Computer { class Equation_Computer; };

namespace CASPT2_ALT{

class CASPT2_ALT {
  
  public:
   
   
    std::shared_ptr<CASPT2::CASPT2> orig_cpt2;
    std::shared_ptr<const SMITH_Info<double>> ref;

    CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in);
    std::shared_ptr<const Dvec> cc_; 
    int  nelea_ ;
    int  neleb_ ;
    int  ncore_ ;
    int  norb_  ;
    int  nstate_;
    std::shared_ptr<const Determinants> det_ ; 

    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> gamma_data_map;
  
    std::shared_ptr<IndexRange> closed_rng; 
    std::shared_ptr<IndexRange> active_rng;  
    std::shared_ptr<IndexRange> virtual_rng;
    std::shared_ptr<IndexRange> free_rng  ; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map ;

    std::vector<std::shared_ptr<MultiTensor_<double>>> T2_all;
    std::vector<std::shared_ptr<MultiTensor_<double>>> lambda_all;
    std::shared_ptr<Tensor_<double>> H_1el_all;
    std::shared_ptr<Tensor_<double>> H_2el_all;// only {occ, virt, occ, virt});

    CASPT2_ALT(std::shared_ptr<const SMITH_Info<double>> ref_alt);
    ~CASPT2_ALT() {};
    
    void test();
    void build_data_map();

    static std::string flip(std::string idx);
    static std::shared_ptr<std::vector<std::string>> ijkl_to_klij(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> ijkl_to_jilk(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> ijkl_to_lkji(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> ijkl_to_ijlk_block(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> ijkl_to_jikl_block(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> ijkl_to_jilk_block(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> bbbb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> bbaa_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> aabb_to_aaaa(std::shared_ptr<std::vector<std::string>> invec) ;
    static std::shared_ptr<std::vector<std::string>> identity(std::shared_ptr<std::vector<std::string>> invec) ;
    static bool NotAllAct(std::shared_ptr<std::vector<std::string>> ranges);
    static bool always_true(std::shared_ptr<std::vector<std::string>> ranges);
    std::vector<std::tuple<std::shared_ptr<std::vector<std::string>>(*)(std::shared_ptr<std::vector<std::string>>),int,int >> set_2el_symmfuncs();

};
}//end CASPT2_ALT namespace

}//end SMITH namespace 
}//end bagel namespace
#endif

