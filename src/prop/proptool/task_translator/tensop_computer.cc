#include <bagel_config.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::Tensor_Sorter;
using namespace bagel::Tensor_Arithmetic; 
using namespace bagel::Tensor_Arithmetic_Utils; 
using namespace WickUtils;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
TensOp_Computer::TensOp_Computer<DataType>::TensOp_Computer( shared_ptr< map< string, shared_ptr<vector<shared_ptr<CtrOp_base>>>>> ACompute_map_in,
                                                             shared_ptr< map< string, shared_ptr<CtrTensorPart_Base>>> CTP_map_in,
                                                             shared_ptr< map< string, shared_ptr<IndexRange>>> range_conversion_map,
                                                             shared_ptr< map< string, shared_ptr<Tensor_<DataType>>>> Data_map_in ):
                                                             ACompute_map(ACompute_map_in), CTP_map(CTP_map_in), Data_map(Data_map_in),
                                                             range_conversion_map_(range_conversion_map), 
                                                             Tensor_Calc(make_shared<Tensor_Arithmetic::Tensor_Arithmetic<DataType>>()){}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
TensOp_Computer::TensOp_Computer<DataType>::Calculate_CTP( AContribInfo& AInfo ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " TensOp_Computer::TensOp_Computer::Calculate_CTP : " <<  AInfo.name_  << endl;

  string A_contrib = AInfo.name_;

  if (ACompute_map->at(A_contrib)->size() == 0 )
    cout << "THIS COMPUTE LIST IS EMPTY" << endl;

  for (shared_ptr<CtrOp_base> ctr_op : *(ACompute_map->at(A_contrib))){ cout << "getting " <<  ctr_op->Tout_name() << endl;
   
    // check if this is an uncontracted multitensor (0,0) && check if the data is in the map
    if( Data_map->find(ctr_op->Tout_name()) == Data_map->end() ) { cout << A_contrib << " not in data_map, must calculate it " << endl;
       shared_ptr<Tensor_<DataType>>  New_Tdata; 
 
      if ( ctr_op->ctr_type()[0] == 'g'){  cout << " : no contraction, fetch this tensor part" << endl; 
        New_Tdata =  get_block_Tensor(ctr_op->Tout_name());
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
        cout << ctr_op->Tout_name() << "->norm() = " << New_Tdata->norm() << endl;
  
      } else if ( ctr_op->ctr_type()[0] == 'd' ){ cout << " : contract different tensors" << endl; 
        New_Tdata = contract_different_tensors( ctr_op->T1name(), ctr_op->T2name(), ctr_op->Tout_name(), make_pair(ctr_op->T1_ctr_rel_pos(), ctr_op->T2_ctr_rel_pos()) );
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
        cout << ctr_op->Tout_name() << "->norm() = " << New_Tdata->norm() << endl;
      
      } else if ( ctr_op->ctr_type()[0] == 's' ) { cout << " : contract on same tensor" <<  endl; 
        New_Tdata = contract_on_same_tensor( ctr_op->T1name(), ctr_op->Tout_name(), ctr_op->ctr_rel_pos() ); 
        Data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
        cout << ctr_op->Tout_name() << "->norm() = " << New_Tdata->norm() << endl;
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
//Returns a tensor T3 with elements T3_{ijkl..} = T1_{ijkl..}/T2_{ijkl..} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> TensOp_Computer::TensOp_Computer<DataType>::divide_tensors(string T1_name, string T2_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " TensOp_Computer::TensOp_Computer::divide_tensors " << endl;
 
   return Tensor_Calc->divide_tensors( find_or_get_CTP_data(T1_name),  find_or_get_CTP_data(T2_name)); 

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//As above, but does division in place, modifying the input;  T1_{ijkl..} = T1_{ijkl..}/T2_{ijkl..} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::divide_tensors_in_place(string T1_name, string T2_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " TensOp_Computer::TensOp_Computer::divide_tensors_in_place" << endl;
 
  Tensor_Calc->divide_tensors_in_place( find_or_get_CTP_data(T1_name),  find_or_get_CTP_data(T2_name)); 

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> TensOp_Computer::TensOp_Computer<DataType>::get_block_Tensor(string Tname){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "TensOp_Computer::TensOp_Computer::get_block_Tensor : " << Tname << endl;
 
   shared_ptr<Tensor_<DataType>> tens; 
   
   if(  Data_map->find(Tname) != Data_map->end()){
     tens = Data_map->at(Tname);

   } else {
     shared_ptr<vector<IndexRange>> id_block = Get_Bagel_IndexRanges( CTP_map->at(Tname)->unc_id_ranges() ) ;

     if(  Data_map->find(Tname.substr(0,1)) != Data_map->end()){
       cout << "initializing uncontracted tensor block " << Tname << " using data from parent tensor \"" << Tname.substr(0,1) << "\"" << endl;
       print_vector( *(CTP_map->at(Tname)->unc_id_ranges()) , "unc_id_ranges" ) ; cout <<endl;  
       tens = get_sub_tensor( Data_map->at(Tname.substr(0,1)), *id_block );

       cout << Tname<< "->norm() = " << tens->norm() << endl; 
 
     } else if ( Tname[0] == 'X' || Tname[0] == 'T' ) {  
       cout << "new tensor block : " << Tname << " is being initialized to 1.0,  for exciations" << endl; 
       tens = make_shared<Tensor_<DataType>>(*id_block);
       tens->allocate();
       Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( tens, 1.0);

     } else {  
       cout << "new tensor block : " << Tname << " is being initialized to zero" << endl; 
       tens = make_shared<Tensor_<DataType>>(*id_block);
       tens->allocate();
       tens->zero();
     }
   }

   cout << "leaving TensOp_Computer::TensOp_Computer::get_block_Tensor : " << Tname << endl;
   return tens;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with ranges specified by unc_ranges, where all values are equal to XX  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::get_uniform_Tensor(shared_ptr<vector<string>> unc_ranges, DataType XX ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "TensOp_Computer::TensOp_Computer::get_uniform_Tensor" << endl;

   shared_ptr<vector<IndexRange>> T_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   return Tensor_Calc->get_uniform_Tensor(T_id_ranges, XX);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_on_same_tensor( std::string T_in_name, std::string T_out_name, pair<int,int> ctr_todo) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "TensOp_Computer::TensOp_Computer::contract_on_same_tensor "; cout.flush();
   cout << ": "  << T_in_name << " over (" << ctr_todo.first << ", " << ctr_todo.second << ") to get " << T_out_name <<  endl;
   return Tensor_Calc->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), ctr_todo  ); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_on_same_tensor( string T_in_name, shared_ptr<vector<int>> ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "TensOp_Computer::TensOp_Computer::contract_on_same_tensor" << ": "  << T_in_name << " over "; print_vector( *ctrs_pos) ; cout << endl;
   return Tensor_Calc->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), *ctrs_pos  ); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                                  pair<int,int> ctr_todo ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_Computer::contract_on_different_tensor" << ": "  << T1_in_name << " and " << T2_in_name << " over T1[" << ctr_todo.first << "] and T2[" << ctr_todo.second << "] to get " << T_out_name <<  endl;
 
  shared_ptr<Tensor_<DataType>> Tens1_in = find_or_get_CTP_data(T1_in_name);
  shared_ptr<Tensor_<DataType>> Tens2_in = find_or_get_CTP_data(T2_in_name);
  return Tensor_Calc->contract_different_tensors(Tens1_in, Tens2_in, ctr_todo );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                              std::pair<std::vector<int>,std::vector<int>> ctrs_todo                  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp_Computer::contract_on_different_tensor" << endl; 

  shared_ptr<Tensor_<DataType>> Tens1_in = find_or_get_CTP_data(T1_in_name);
  shared_ptr<Tensor_<DataType>> Tens2_in = find_or_get_CTP_data(T2_in_name);

  return Tensor_Calc->contract_different_tensors(Tens1_in, Tens2_in, ctrs_todo );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::direct_product_tensors( std::vector<std::string>& Tensor_names ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp_Computer::direct_product_tensors" << endl; 

  print_vector(Tensor_names, "Tensor_names"); cout <<endl; 

  shared_ptr<Tensor_<DataType>> Tens_prod = find_or_get_CTP_data(Tensor_names[0]);
  shared_ptr<Tensor_<DataType>> Tens_intermediate;
  for ( int rr = 1 ; rr != Tensor_names.size() ; rr++) {
    shared_ptr<Tensor_<DataType>> Tens_next = find_or_get_CTP_data(Tensor_names[rr]);
    Tens_intermediate = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::direct_tensor_product( Tens_prod, Tens_next ); 
    Tens_prod = Tens_intermediate;
  }           
 
  return Tens_prod;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::reorder_block_Tensor(string T_in_name, shared_ptr<vector<int>> new_order){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "TensOp_Computer::TensOp_Computer::reorder_block_Tensor "; cout.flush();
  cout << " : " << T_in_name ; cout.flush();
  cout <<  " New_order = [ "; cout.flush();  for (int pos : *new_order ) { cout << pos << " " ; cout.flush(); } cout << "] " << endl;
  
  auto Data_map_loc = Data_map->find(T_in_name); 
  shared_ptr<Tensor_<DataType>> T_part; 
  if( Data_map_loc == Data_map->end() ){
    T_part =  get_block_Tensor(T_in_name); 
    Data_map->emplace(T_in_name , T_part );
  } else { 
    T_part = Data_map_loc->second; 
  }

  return Tensor_Calc->reorder_block_Tensor( T_part , new_order );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
pair<int,int>
TensOp_Computer::TensOp_Computer<DataType>::relativize_ctr_positions( pair <int,int> ctr_todo, 
                                                            shared_ptr<CtrTensorPart_Base> CTP1,
                                                            shared_ptr<CtrTensorPart_Base> CTP2 ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "TensOp_Computer::TensOp_Computer::relativize_ctr_positions" << endl;
   pair<int,int> rel_ctr;

   int T1_orig_size = CTP1->full_id_ranges()->size(); 
   int T2_orig_size = CTP2->full_id_ranges()->size(); 

   if (ctr_todo.first >= T1_orig_size ){ 
     rel_ctr = make_pair(ctr_todo.second, ctr_todo.first-T2_orig_size);
   } else {
     rel_ctr = make_pair(ctr_todo.first, ctr_todo.second-T1_orig_size);
   }

  for ( int ii = 0 ; ii != CTP1->unc_pos()->size() ; ii++ )
     if ( CTP1->unc_pos(ii) == rel_ctr.first){
       rel_ctr.first = ii;
       break;
     }

  for ( int ii = 0 ; ii != CTP2->unc_pos()->size() ; ii++ )
     if (CTP2->unc_pos(ii) == rel_ctr.second){ 
       rel_ctr.second = ii;
       break;
     }

 return rel_ctr;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO This should be replaced; ultimately, each of the blocks should be it's own tensor; this duplicates everything.
//     Which blocks exist should be determined after the equation editor has run, so this should ultimately link to
//     allocation and initialization routines, etc..
//     As an intermediate step could perhaps erase the newly created CTP block as soon as it's extracted from the tensor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::find_or_get_CTP_data(string CTP_name){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 cout << "TensOp_Computer::TensOp_Computer::find_or_get_CTP_data  : " <<  CTP_name <<  endl;

  shared_ptr<Tensor_<DataType>> CTP_data;
  auto Data_loc =  Data_map->find(CTP_name); 
  if ( Data_loc == Data_map->end() ){
    CTP_data = get_block_Tensor(CTP_name);   
    Data_map->emplace(CTP_name, CTP_data);   
  } else {
    CTP_data = Data_loc->second;
  }
  
  return CTP_data;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<vector<shared_ptr<const IndexRange>>>
TensOp_Computer::TensOp_Computer<DataType>::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 1arg" << endl;

  shared_ptr<vector<shared_ptr<const IndexRange>>> ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map_->at(rng)));

  return ranges_Bagel;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<vector<shared_ptr<const IndexRange>>>
TensOp_Computer::TensOp_Computer<DataType>::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str, shared_ptr<vector<int>> unc_pos){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_Computer::Get_Bagel_const_IndexRanges 2arg" ; print_vector(*ranges_str, "ranges_str" ) ;
cout <<  "  "; cout.flush();  print_vector(*unc_pos, "unc_pos" ) ; cout << endl;

  vector<shared_ptr<const IndexRange>>  ranges_Bagel(unc_pos->size());
  for ( int ii = 0 ; ii != unc_pos->size() ; ii++) 
    ranges_Bagel[ii]=(make_shared<const IndexRange>(*range_conversion_map_->at(ranges_str->at(unc_pos->at(ii)))));

  return make_shared<vector<shared_ptr<const IndexRange>>>(ranges_Bagel);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<vector<IndexRange>>
TensOp_Computer::TensOp_Computer<DataType>::Get_Bagel_IndexRanges( shared_ptr<vector<string>> ranges_str ){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "TensOp_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(*ranges_str, "ranges_str" ) ; cout << endl;

  shared_ptr<vector<IndexRange>> ranges_Bagel = make_shared<vector<IndexRange>>(ranges_str->size());
  for ( int ii = 0 ; ii != ranges_str->size(); ii++)
    ranges_Bagel->at(ii) = *range_conversion_map_->at(ranges_str->at(ii));

  return ranges_Bagel;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp_Computer::TensOp_Computer<double>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
