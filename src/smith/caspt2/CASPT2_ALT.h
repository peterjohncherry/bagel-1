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
#include <src/smith/wicktool/expression_info.h>
#include <src/smith/caspt2/CASPT2.h>
#include <src/ci/fci/fci.h>
#include <src/smith/wicktool/equation_computer.h>
#include <src/smith/wicktool/gamma_computer.h>
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

    std::shared_ptr<const Dvec> cc_; 
    int  ncore ;
    int  nact  ;
    int  nvirt ;
    int  nclosed ;
    int  nstate;
    std::shared_ptr<const Determinants> det_ ; 

    std::shared_ptr<std::map< std::string, std::shared_ptr<CtrTensorPart<double>>>> CTP_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor_<double>>>> Gamma_data_map;
    std::shared_ptr<std::map< std::string, std::shared_ptr<Equation<double>>>> Expr_Info_map;
    std::shared_ptr<std::map< std::string, double>> scalar_results_map;
  
    std::shared_ptr<IndexRange> closed_rng; 
    std::shared_ptr<IndexRange> active_rng;  
    std::shared_ptr<IndexRange> virtual_rng;
    std::shared_ptr<IndexRange> free_rng  ; 
    std::shared_ptr<IndexRange> not_closed_rng ; 
    std::shared_ptr<IndexRange> not_active_rng  ;
    std::shared_ptr<IndexRange> not_virtual_rng ;

    std::shared_ptr<StatesInfo<double>> TargetsInfo;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map ;
     
    std::vector<std::shared_ptr<MultiTensor_<double>>> T2_all;
    std::vector<std::shared_ptr<MultiTensor_<double>>> lambda_all;
    std::shared_ptr<Tensor_<double>> F_1el_all;
    std::shared_ptr<Tensor_<double>> H_1el_all;
    std::shared_ptr<Tensor_<double>> H_2el_all;// only {occ, virt, occ, virt});

    std::shared_ptr<Expression_Info<double>> Expr_Info;

    CASPT2_ALT(const CASPT2::CASPT2& orig_cpt2_in);
    ~CASPT2_ALT() {};
    
    void solve();
    void Construct_Tensor_Ops();
    void build_data_map();
    void Build_Compute_Lists();
    void Execute_Compute_List(std::string expression_name );
    
    void set_target_info( std::shared_ptr<std::vector<int>> states_of_interest) ;

    struct compare_string_length {
      bool operator()(const std::string& first, const std::string& second) {
          return first.size() > second.size();
      }
};


};
}//end CASPT2_ALT namespace

}//end SMITH namespace 
}//end bagel namespace
#endif

