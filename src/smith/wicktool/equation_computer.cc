#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_computer.h>
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;
using namespace bagel::SMITH::Tensor_Arithmetic; 
using namespace bagel::SMITH::Tensor_Arithmetic_Utils; 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Equation_Computer::Equation_Computer::Equation_Computer( std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<double>> eqn_info_in,
                                                         std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                                                         std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> Data_map_in          ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Equation_Computer::Equation_Computer::Equation_Computer" << endl;
  eqn_info =  eqn_info_in;
  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  cc_ = ref->ciwfn()->civectors();
  maxtile = ref->maxtile();
  
  CTP_map  = eqn_info->CTP_map;
  GammaMap = eqn_info->GammaMap;
  
  Data_map  = Data_map_in;

  cimaxblock = 100; //figure out what is best, maxtile is 10000, so this is chosen to have one index block. Must be consistent if contraction routines are to work...

  range_conversion_map = range_conversion_map_in;

  Tensor_Calc = make_shared<Tensor_Arithmetic::Tensor_Arithmetic<double>>();

}  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::Calculate_CTP(std::string A_contrib ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " Equation_Computer::Equation_Computer::Calculate_CTP(std::string A_contrib_name )" << endl;
  if (eqn_info->ACompute_map->at(A_contrib)->size() == 0 )
    cout << "THIS COMPUTE LIST IS EMPTY" << endl;

  for (shared_ptr<CtrOp_base> ctr_op : *(eqn_info->ACompute_map->at(A_contrib))){
 
    // check if this is an uncontracted multitensor (0,0) && check if the data is in the map
    if( Data_map->find(ctr_op->Tout_name()) == Data_map->end() ) {
       shared_ptr<Tensor_<double>>  New_Tdata; 
 
      if ( ctr_op->ctr_type()[0] == 'g'){  cout << " : no contraction, fetch this tensor part" << endl; 
        New_Tdata =  get_block_Tensor(ctr_op->Tout_name());
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
  
      } else if ( ctr_op->ctr_type()[0] == 'd' ){ cout << " : contract different tensors" << endl; 
        New_Tdata = contract_different_tensors( ctr_op->T1name(), ctr_op->T2name(), ctr_op->Tout_name(), make_pair(ctr_op->T1_ctr_rel_pos(), ctr_op->T2_ctr_rel_pos()));
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
      
      } else if ( ctr_op->ctr_type()[0] == 's' ) { cout << " : contract on same tensor" <<  endl; 
        New_Tdata = contract_on_same_tensor( ctr_op->T1name(), ctr_op->Tout_name(), ctr_op->ctr_rel_pos() ); 
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
      } else { 
        throw std::runtime_error(" unknown contraction type : " + ctr_op->ctr_type() ) ;
      }
    }
    cout << "A_contrib.first = " << A_contrib << endl;
    cout << "CTPout_name =  " << ctr_op->Tout_name() << endl;
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_block_Tensor(string Tname){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_block_Tensor : " << Tname << endl;
  
   if(  Data_map->find(Tname) != Data_map->end())
     return Data_map->at(Tname);

   shared_ptr<vector<string>> unc_ranges = CTP_map->at(Tname)->unc_id_ranges;  

   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   shared_ptr<vector<int>> range_lengths = get_range_lengths( Bagel_id_ranges ) ;

   shared_ptr<Tensor_<double>> fulltens = Data_map->at(Tname.substr(0,1));
   shared_ptr<Tensor_<double>> block_tensor = make_shared<Tensor_<double>>(*Bagel_id_ranges);
   block_tensor->allocate();
   block_tensor->zero();

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges->size(),0);  

   do {
     cout << Tname << " block pos =  [ " ;    for (int block_num : *block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
     
     vector<Index> T_id_blocks =  *(get_rng_blocks( block_pos, *Bagel_id_ranges )); 

     unique_ptr<double[]> T_block_data = fulltens->get_block(T_id_blocks);

     block_tensor->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));

   return block_tensor;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with ranges specified by unc_ranges, where all values are equal to XX  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_uniform_Tensor(shared_ptr<vector<string>> unc_ranges, double XX ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_uniform_Tensor" << endl;

   shared_ptr<vector<IndexRange>> T_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   return Tensor_Calc->get_uniform_Tensor(T_id_ranges, XX);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_on_same_tensor( std::string T_in_name, std::string T_out_name, pair<int,int> ctr_todo) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::contract_on_same_tensor "; cout.flush();
   cout << ": "  << T_in_name << " over (" << ctr_todo.first << ", " << ctr_todo.second << ") to get " << T_out_name <<  endl;

   return Tensor_Calc->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), ctr_todo  ); 

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::contract_on_same_tensor( string T_in_name, shared_ptr<vector<int>> ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::contract_on_same_tensor" << endl;
   cout << ": "  << T_in_name << " over ( "; cout.flush();
   for ( int pos : *ctrs_pos){ 
     cout << pos << " " ; cout.flush();
   }
   cout << " ) " << endl;
   return Tensor_Calc->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), *ctrs_pos  ); 

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                                  pair<int,int> ctr_todo ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::contract_on_different_tensor" <<endl; 
cout << ": "  << T1_in_name << " and " << T2_in_name << " over T1[" << ctr_todo.first << "] and T2[" << ctr_todo.second << "] to get " << T_out_name <<  endl;

  shared_ptr<Tensor_<double>> Tens1_in = find_or_get_CTP_data(T1_in_name);
  shared_ptr<Tensor_<double>> Tens2_in = find_or_get_CTP_data(T2_in_name);

  shared_ptr<Tensor_<double>> Tens_out = Tensor_Calc->contract_different_tensors(Tens1_in, Tens2_in, ctr_todo );
  cout << "CTP1_data->rms()     = " + T1_in_name + "->rms() =" << Tens1_in->rms() << endl;
  cout << "CTP2_data->rms()     = " + T2_in_name + "->rms() =" << Tens2_in->rms() << endl;
  cout << "CTP_data_out->rms()  = " + T_out_name + "->rms() =" << Tens_out->rms() << endl;
 
  return Tens_out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::reorder_block_Tensor(string T_in_name, shared_ptr<vector<int>> new_order){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::reorder_block_Tensor "; cout.flush();
   cout << " : " << T_in_name ; cout.flush();
   cout <<  " New_order = [ "; cout.flush();  for (int pos : *new_order ) { cout << pos << " " ; cout.flush(); } cout << "] " << endl;
  
   auto Data_map_loc = Data_map->find(T_in_name); 
   if(  Data_map_loc == Data_map->end()){
     throw std::runtime_error(" don't have tensor data for " +T_in_name+ " yet.... Equation_Computer::reorder_block_Tensor " ) ;
   } else { 
     return Tensor_Calc->reorder_block_Tensor( Data_map_loc->second , new_order);
   }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<int,int>
Equation_Computer::Equation_Computer::relativize_ctr_positions( pair <int,int> ctr_todo, 
                                                                shared_ptr<CtrTensorPart<double>> CTP1,
                                                                shared_ptr<CtrTensorPart<double>> CTP2 ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Equation_Computer::Equation_Computer::find_or_get_CTP_data" << endl;
   pair<int,int> rel_ctr;

   int T1_orig_size = CTP1->full_id_ranges->size(); 
   int T2_orig_size = CTP2->full_id_ranges->size(); 

   if (ctr_todo.first >= T1_orig_size ){ 
     rel_ctr = make_pair(ctr_todo.second, ctr_todo.first-T2_orig_size);
   } else {
     rel_ctr = make_pair(ctr_todo.first, ctr_todo.second-T1_orig_size);
   }

  for ( int ii = 0 ; ii != CTP1->unc_pos->size() ; ii++ )
     if ( CTP1->unc_pos->at(ii) == rel_ctr.first){
       rel_ctr.first = ii;
       break;
     }

  for ( int ii = 0 ; ii != CTP2->unc_pos->size() ; ii++ )
     if (CTP2->unc_pos->at(ii) == rel_ctr.second){ 
       rel_ctr.second = ii;
       break;
     }

 return rel_ctr;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::find_or_get_CTP_data(string CTP_name){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "Equation_Computer::Equation_Computer::find_or_get_CTP_data  : " <<  CTP_name <<  endl;

  shared_ptr<Tensor_<double>> CTP_data;
  auto Data_loc =  Data_map->find(CTP_name); 
  if ( Data_loc == Data_map->end() ){
    CTP_data = get_block_Tensor(CTP_name);   
  } else {
    CTP_data = Data_loc->second;
  }
  
  return CTP_data;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 1arg" << endl;

  shared_ptr<vector<shared_ptr<const IndexRange>>> ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map->at(rng)));

  return ranges_Bagel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str, shared_ptr<vector<int>> unc_pos){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::Get_Bagel_const_IndexRanges 2arg" ; print_vector(*ranges_str, "ranges_str" ) ;
cout <<  "  "; cout.flush();  print_vector(*unc_pos, "unc_pos" ) ; cout << endl;

  vector<shared_ptr<const IndexRange>>  ranges_Bagel(unc_pos->size());
  for ( int ii = 0 ; ii != unc_pos->size() ; ii++) 
    ranges_Bagel[ii]=(make_shared<const IndexRange>(*range_conversion_map->at(ranges_str->at(unc_pos->at(ii)))));

  return make_shared<vector<shared_ptr<const IndexRange>>>(ranges_Bagel);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> Equation_Computer::Equation_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(*ranges_str, "ranges_str" ) ; cout << endl;

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}
#endif
