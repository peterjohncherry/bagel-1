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

  civectors = orig_cpt2_in.info_->ciwfn()->civectors();

  T2_all     = orig_cpt2_in.t2all_;  
  lambda_all = orig_cpt2_in.lall_;   
  F_1el_all  = orig_cpt2_in.f1_;     cout << "F_1el_all->rms()" <<  F_1el_all->rms() << endl; 
  H_1el_all  = orig_cpt2_in.h1_;     cout << "H_1el_all->rms()" <<  H_1el_all->rms() << endl; 
  H_2el_all  = orig_cpt2_in.H_2el_;  cout << "H_2el_all->rms()" <<  H_2el_all->rms() << endl; 
 
  ncore   = orig_cpt2_in.info_->ncore();
  nclosed = orig_cpt2_in.info_->nclosed();
  nact    = orig_cpt2_in.info_->nact();
  nvirt   = orig_cpt2_in.info_->nvirt();
  maxtile = orig_cpt2_in.info_->maxtile();
  cimaxtile = 100000;

  TensOp_data_map    = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  Gamma_data_map     = make_shared<map<string, shared_ptr<Tensor_<double>>>>();
  scalar_results_map = make_shared<map<string, double>>();
  range_conversion_map = make_shared<map<string, shared_ptr<IndexRange>>>();

  CIvec_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    
  Sigma_data_map   = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();    

  shared_ptr<vector<int>> states_of_interest = make_shared<vector<int>>( vector<int> { 0 } );
  set_target_info(states_of_interest) ;
  set_range_info(states_of_interest);
   
  Sys_Info = make_shared<System_Info<double>>(TargetsInfo, true);
  Expression_map = Sys_Info->expression_map;
  Expression_Machine = make_shared<Expression_Computer::Expression_Computer<double>>(  civectors, Expression_map, range_conversion_map, TensOp_data_map, 
                                                                                       Gamma_data_map, Sigma_data_map, CIvec_data_map );

  ////////////////// FOR TESTING ///////////////////////////////////////  
  shared_ptr<CASPT2::CASPT2> orig_cpt2 = make_shared<CASPT2::CASPT2>(orig_cpt2_in);
  orig_cpt2->set_rdm(0,0);
  Smith_rdm1 = orig_cpt2->rdm1_; TensOp_data_map->emplace( "Smith_rdm1", Smith_rdm1 );
  Smith_rdm2 = orig_cpt2->rdm2_; TensOp_data_map->emplace( "Smith_rdm2", Smith_rdm2 );
  Smith_rdm3 = orig_cpt2->rdm3_; TensOp_data_map->emplace( "Smith_rdm3", Smith_rdm3 );
  Smith_rdm4 = orig_cpt2->rdm4_; TensOp_data_map->emplace( "Smith_rdm4", Smith_rdm4 );

  shared_ptr<vector<int>> gorder = make_shared<vector<int>>( vector<int> {  1, 2, 0, 3} ); 
  shared_ptr<Tensor_<double>> rdm2_reord = 
  Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor( Smith_rdm2, gorder );
 
  vector<Index> id_block_rdm2 = { active_rng->range(0), active_rng->range(0), active_rng->range(0),active_rng->range(0)}; 
  vector<Index> id_block_rdm1 = { active_rng->range(0), active_rng->range(0) }; 

  unique_ptr<double[]> rdm2_data = rdm2_reord->get_block( id_block_rdm2 );
  unique_ptr<double[]> rdm1_data = Smith_rdm1->get_block( id_block_rdm1 );
 
}
//////////////////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::set_range_info(shared_ptr<vector<int>> states_of_interest ) {
//////////////////////////////////////////////////////////////////////////////////////////////
  cout << "CASPT2_ALT::CASPT2_ALT::set_range_info " << endl;

  cout << maxtile << endl;
  int maxtile_buff = maxtile; 
  closed_rng  =  make_shared<IndexRange>(IndexRange(nclosed-ncore, maxtile, 0, ncore));
  active_rng  =  make_shared<IndexRange>(IndexRange(nact, min(10,maxtile_buff), closed_rng->nblock(), ncore + closed_rng->size()));
  virtual_rng =  make_shared<IndexRange>(IndexRange(nvirt, maxtile, closed_rng->nblock()+ active_rng->nblock(), ncore+closed_rng->size()+active_rng->size()));
  free_rng    = make_shared<IndexRange>(*closed_rng);
  free_rng->merge(*active_rng);
  free_rng->merge(*virtual_rng);

  not_closed_rng  =  make_shared<IndexRange>(*active_rng); not_closed_rng->merge(*virtual_rng);
  not_active_rng  =  make_shared<IndexRange>(*closed_rng); not_active_rng->merge(*virtual_rng);
  not_virtual_rng =  make_shared<IndexRange>(*closed_rng); not_virtual_rng->merge(*active_rng);

  PT2_ranges_ = { *not_closed_rng, *not_closed_rng, *not_virtual_rng, *not_virtual_rng};
  PT2_ranges_herm_conj_ = {  *not_virtual_rng, *not_virtual_rng, *not_closed_rng, *not_closed_rng};

  range_conversion_map->emplace("cor", closed_rng);//change the naming of the ranges from cor to clo... 
  range_conversion_map->emplace("act", active_rng);
  range_conversion_map->emplace("A", active_rng); //fudge, but leave for now
  range_conversion_map->emplace("vir", virtual_rng);
  range_conversion_map->emplace("free", free_rng);

  range_conversion_map->emplace("notcor", not_closed_rng);
  range_conversion_map->emplace("notact", not_active_rng);
  range_conversion_map->emplace("notvir", not_virtual_rng); 
 
  for ( int ii : *states_of_interest ) {

    range_conversion_map->emplace( get_civec_name( ii , civectors->data(ii)->det()->norb(), civectors->data(ii)->det()->nelea(), civectors->data(ii)->det()->neleb()),
                                                   make_shared<IndexRange>(civectors->data(ii)->det()->size(), cimaxtile ));  
    
    shared_ptr<IndexRange>  ci_index_ranges =  make_shared<IndexRange>(civectors->data(ii)->det()->size(), cimaxtile );
    cout << "cirngs = [ ";  for (auto irng : ci_index_ranges->range()) { cout << irng.size()  << " "; }; cout << "] " << endl;

  }    

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
  TargetsInfo = make_shared<StatesInfo<double>> ( *states_of_interest ) ;
  
  for ( int state_num : *states_of_interest ) 
     TargetsInfo->add_state( civectors->data(state_num)->det()->nelea(), civectors->data(state_num)->det()->neleb(),
                             civectors->data(state_num)->det()->norb(), state_num ) ;
  
  return;
}

////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::solve() {
////////////////////////////////////////////////////////////////////
cout <<  " CASPT2_ALT::CASPT2_ALT::solve() " << endl;

 {
   vector<string> op_list = { "L", "Q" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LQ" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }

 {
   vector<string> op_list = { "M", "N" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "MN" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }

 {
   vector<string> op_list = { "L", "R" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LR" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }

 {
   vector<string> op_list = { "L", "S" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LS" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }

 {
   vector<string> op_list = { "L", "T" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LT" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }


 {
   vector<string> op_list = { "L", "U" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LU" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }


 {
   vector<string> op_list = { "L", "V" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LV" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }


 {
   vector<string> op_list = { "L", "W" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LW" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }


 {
   vector<string> op_list = { "L", "Y" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( op_list, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LY" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
     }
   }
 }

 {
   vector<string> lz = { "L", "Z" };
   vector< pair<vector<string>,double> > BK_info_list( 1, make_pair( lz, 1.0 ) );
   
   double factor = 0.0;
   // Building all necessary expressions 
   int  num_states = 1; 
   vector<vector<Term_Info<double>>> Term_info_list( num_states*num_states );
   for ( int ii = 0 ; ii != num_states; ii++) {
     for ( int jj = 0 ; jj != num_states; jj++) {
   
       for ( pair<vector<string>,double> BK_info : BK_info_list ) {
         Term_info_list[ii*num_states+jj].push_back(Term_Info<double>( BK_info.first, TargetsInfo->name(ii), TargetsInfo->name(jj), BK_info.second , "LZ" ));
         for ( string Op_name : BK_info.first )  
           Set_Tensor_Ops_Data( Op_name, TargetsInfo->name(ii), TargetsInfo->name(jj) ); 
       }
   
       string expression_name = Sys_Info->Build_Expression( Term_info_list[ii*num_states+jj] );
   
       Expression_Machine->Evaluate_Expression( expression_name );
  
      }
    }
  }


 return;
} 
/////////////////////////////////////////////////////////////////////////////////
//Build the operators here. 
/////////////////////////////////////////////////////////////////////////////////
void CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data( string op_name, string bra_name , string ket_name )  { 
/////////////////////////////////////////////////////////////////////////////////
cout << "CASPT2_ALT::CASPT2_ALT::Set_Tensor_Ops_Data() " << endl;

  vector<IndexRange> ffff = { *free_rng, *free_rng, *free_rng, *free_rng};

  vector<IndexRange> aaaa = { *active_rng,  *active_rng,  *active_rng,      *active_rng};
  vector<IndexRange> vvaa = { *virtual_rng, *virtual_rng, *active_rng,      *active_rng} ;

  vector<IndexRange> vvcc = { *virtual_rng, *virtual_rng, *closed_rng,      *closed_rng} ;
  vector<IndexRange> vvac = { *virtual_rng, *virtual_rng, *active_rng,      *closed_rng} ;
  vector<IndexRange> vvca = { *virtual_rng, *virtual_rng, *closed_rng,      *active_rng} ;

  vector<IndexRange> vvoo = { *virtual_rng,     *virtual_rng,     *not_virtual_rng, *not_virtual_rng } ;
  vector<IndexRange> vovo = { *virtual_rng,     *not_virtual_rng, *virtual_rng,     *not_virtual_rng } ;
  vector<IndexRange> ovov = { *not_virtual_rng, *virtual_rng,     *not_virtual_rng, *virtual_rng     } ;

  vector<IndexRange> vvoc = { *virtual_rng, *virtual_rng, *not_virtual_rng, *closed_rng} ;
  vector<IndexRange> vvco = { *virtual_rng, *virtual_rng, *closed_rng,      *not_virtual_rng} ;
  vector<IndexRange> vvoa = { *virtual_rng, *virtual_rng, *not_virtual_rng, *active_rng} ;
  vector<IndexRange> vvao = { *virtual_rng, *virtual_rng, *active_rng,      *not_virtual_rng} ;

  //Setting_data for TensOps
  if ( op_name  == "P" ) { 

    shared_ptr<vector<IndexRange>> H2el_ranges = make_shared<vector<IndexRange>>( H_2el_all->indexrange() );
    shared_ptr<Tensor_<double>> P_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( H2el_ranges);
    TensOp_data_map->emplace( "P", P_Tens );

  } else if ( op_name  == "Smith_rdms" ) { 
   
    TensOp_data_map->emplace("Smith_rdm1", Smith_rdm1 );
    TensOp_data_map->emplace("Smith_rdm2", Smith_rdm2 );

  } else if ( op_name  == "F" ) {
 
    shared_ptr<Tensor_<double>> F_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( make_shared<vector<IndexRange>>(ffff), 1.0 );
    TensOp_data_map->emplace( "F", F_Tens );


  } else if ( op_name  == "G" ) {
 
    shared_ptr<vector<IndexRange>> act2_ranges = make_shared<vector<IndexRange>>( vector<IndexRange> { *not_closed_rng, *not_virtual_rng } );
    shared_ptr<Tensor_<double>> G_Tens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( act2_ranges, 2.0 );
    TensOp_data_map->emplace( "G", G_Tens );

  } else if ( op_name  == "H" ) { 
    
   shared_ptr<Tensor_<double>> HTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor( make_shared<vector<IndexRange>>(ffff), 0.0 );
   Tensor_Arithmetic::Tensor_Arithmetic<double>::set_tensor_elems( HTens , vvcc,  0.5  );
   TensOp_data_map->emplace( "H" , HTens );

  cout << "HTens->norm() = "<< HTens->norm() << endl;

  } else if ( op_name  == "h" ) { 

    TensOp_data_map->emplace("h" , H_1el_all);

  } else if ( op_name  == "L" ) { 

    shared_ptr<Tensor_<double>> LTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(make_shared<vector<IndexRange>>(vvoo), 1.0 );
    TensOp_data_map->emplace("L" , LTens);

  } else if ( op_name  == "M" ) { 

    shared_ptr<Tensor_<double>> MTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_uniform_Tensor(make_shared<vector<IndexRange>>(vovo), 1.0 );
    TensOp_data_map->emplace("M" , MTens);

  } else if ( op_name  == "N" ) { 

    shared_ptr<Tensor_<double>> NTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major(make_shared<vector<IndexRange>>(vovo));
    TensOp_data_map->emplace("N" , NTens);

  } else if ( op_name  == "Q" ) { 

    shared_ptr<Tensor_<double>> QTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major(make_shared<vector<IndexRange>>(vvoo) );
    TensOp_data_map->emplace("Q" , QTens);

  } else if ( op_name  == "R" ) { 

    shared_ptr<Tensor_<double>> RTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvaa) );
    TensOp_data_map->emplace( "R" , RTens);

  } else if ( op_name  == "S" ) { 

    shared_ptr<Tensor_<double>> TTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvao) );
    shared_ptr<vector<int>> vvao_to_vvoa = make_shared<vector<int>>( vector<int> { 0, 1, 3, 2} ) ;  
    shared_ptr<Tensor_<double>> STens = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor(TTens, vvao_to_vvoa);
    TensOp_data_map->emplace( "S", STens );

  } else if ( op_name  == "T" ) { 

    shared_ptr<Tensor_<double>> TTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvao) );
    TensOp_data_map->emplace("T", TTens);

  } else if ( op_name  == "U" ) { 

    shared_ptr<Tensor_<double>> UTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvco) );
    TensOp_data_map->emplace("U", UTens);

  } else if ( op_name  == "V" ) { 

    shared_ptr<Tensor_<double>> UTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvco) );
    shared_ptr<vector<int>> vvco_to_vvoc = make_shared<vector<int>>( vector<int> { 0, 1, 3, 2} ) ;  
    shared_ptr<Tensor_<double>> VTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor(UTens, vvco_to_vvoc);
    TensOp_data_map->emplace("V", VTens);

  } else if ( op_name  == "W" ) { 

    shared_ptr<Tensor_<double>> ZTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major(  make_shared<vector<IndexRange>>(vvac) );
    shared_ptr<vector<int>> vvac_to_vvca = make_shared<vector<int>>( vector<int> { 0, 1, 3, 2} ) ;  
    shared_ptr<Tensor_<double>> WTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::reorder_block_Tensor(ZTens, vvac_to_vvca);
    TensOp_data_map->emplace("W", WTens);

  } else if ( op_name  == "X" ) { 

    shared_ptr<Tensor_<double>> XTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvaa) );
    TensOp_data_map->emplace( "X" , XTens);
 
  } else if ( op_name  == "Y" ) {
 
    shared_ptr<Tensor_<double>> YTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvcc) );
    TensOp_data_map->emplace( "Y", YTens );
   
  } else if ( op_name  == "Z" ) {
 
    shared_ptr<Tensor_<double>> ZTens = Tensor_Arithmetic::Tensor_Arithmetic<double>::get_test_Tensor_column_major( make_shared<vector<IndexRange>>(vvac) );
    TensOp_data_map->emplace( "Z", ZTens );
  }

  return;
}
#endif
