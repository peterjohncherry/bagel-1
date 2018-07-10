#include <bagel_config.h>
#include <src/util/f77.h>
#include <src/prop/proptool/task_translator/tensop_computer.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic_utils.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::Tensor_Arithmetic; 
using namespace bagel::Tensor_Arithmetic_Utils; 
using namespace WickUtils;

#define __DEBUG_TENSOP_COMPUTER
#define __DEBUG_TENSOP_COMPUTER_X
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void
TensOp_Computer::TensOp_Computer<DataType>::Calculate_CTP( AContribInfo_Base& AInfo ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER_X
cout << " TensOp_Computer::TensOp_Computer::Calculate_CTP : "; cout.flush(); cout <<  AInfo.name_ << endl;
print_compute_list( AInfo.name_ );
#endif ////////////////////////////////////////////////////////////////////////////////////////////////////////
  string A_contrib = AInfo.name_;

  if (ACompute_map_->at(A_contrib)->size() == 0 )
    throw logic_error( "ACompute list is empty!! Aborting!!" );

  for (shared_ptr<CtrOp_base> ctr_op : *(ACompute_map_->at(A_contrib))){
   
    //TODO don't check ctr_type(), check if dynamic_pointer_cast<X>( ctr_op ) is null (where X is the type of contraction you're casting to).
    // check if this is an uncontracted multitensor (0,0) && check if the data is in the map
    if( tensop_data_map_->find(ctr_op->Tout_name()) == tensop_data_map_->end() ) {
       cout << A_contrib << " not in data_map, must calculate it " << endl;
       shared_ptr<Tensor_<DataType>>  New_Tdata; 
 
      if ( ctr_op->ctr_type()[0] == 'g'){ 
        New_Tdata =  get_block_Tensor(ctr_op->Tout_name());
        tensop_data_map_->emplace(ctr_op->Tout_name(), New_Tdata); 
  
      } else if ( ctr_op->ctr_type()[0] == 'd' ){
        New_Tdata = contract_different_tensors( ctr_op->T1name(), ctr_op->T2name(), ctr_op->Tout_name(), make_pair(ctr_op->T1_ctr_rel_pos(), ctr_op->T2_ctr_rel_pos()) );
        tensop_data_map_->emplace(ctr_op->Tout_name(), New_Tdata); 
      
      } else if ( ctr_op->ctr_type()[0] == 's' ) {
        New_Tdata = contract_on_same_tensor( ctr_op->T1name(), ctr_op->Tout_name(), ctr_op->ctr_rel_pos() ); 
        tensop_data_map_->emplace(ctr_op->Tout_name(), New_Tdata); 

      } else if ( ctr_op->ctr_type()[0] == 'r' ) {
        New_Tdata = reorder_block_Tensor( ctr_op->T1name(), *(ctr_op->new_order()) ); 
        tensop_data_map_->emplace(ctr_op->Tout_name(), New_Tdata); 

      } else if ( ctr_op->ctr_type()[0] == 'c' ) {
        New_Tdata = direct_product_tensors( *(ctr_op->tensor_list()) );
        tensop_data_map_->emplace(ctr_op->Tout_name(), New_Tdata);

      } else { 
        throw std::runtime_error(" unknown contraction type : " + ctr_op->ctr_type() ) ;
      }

      #ifdef __DEBUG_TENSOP_COMPUTER
         cout << ctr_op->Tout_name() << "->norm() = " << New_Tdata->norm() << endl;
      #endif
    }
    
   #ifdef __DEBUG_TENSOP_COMPUTER
      cout << "A_contrib.first = " << A_contrib; cout.flush(); cout << "  CTPout_name =  " << ctr_op->Tout_name() << endl;
   #endif
  }
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::add_acontrib_to_target( string target_name, const shared_ptr<AContribInfo_Base>& a_contrib ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::add_acontrib_to_target " << endl;
a_contrib->print_info(); cout << endl;
cout << "====================================================================================================" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sum_different_orderings( target_name , a_contrib->name(), a_contrib->factors(), a_contrib->id_orders() );

  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Adds summand name factor to target name 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::sum_different_orderings( string target_name, string summand_name, 
                                                                          vector<pair<double,double>> summand_factors,
                                                                          vector<vector<int>> summand_reorderings  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << endl <<  "TensOp_Computer::TensOp_Computer::sum_different_orderings " << endl;
cout << "target_name  = " << target_name ; cout.flush(); cout << "  summand_name = " << summand_name << endl;
cout << "target->norm() +  summand_orig_order->norm() : "; cout.flush();
cout <<  find_or_get_CTP_data( target_name )->norm(); cout.flush(); cout << " + " << find_or_get_CTP_data( summand_name )->norm(); cout.flush(); 
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> target  = find_or_get_CTP_data( target_name );
  shared_ptr<Tensor_<DataType>> summand_orig_order = find_or_get_CTP_data( summand_name );

  //TODO remove this, and decide properly how you are going to deal with the real and complex factors
  vector<DataType> summand_factors_re( summand_factors.size() );
  vector<DataType> summand_factors_im( summand_factors.size() );
  typename vector<DataType>::iterator sfr_it = summand_factors_re.begin(); 
  for ( vector<pair<double,double>>::iterator sr_it = summand_factors.begin() ; sr_it != summand_factors.end() ; sr_it++ , sfr_it++ )
    *sfr_it = (DataType)(sr_it->first); 

  Tensor_Calc_->add_list_of_reordered_tensors( target, summand_orig_order, summand_reorderings, summand_factors_re );
 
  tensop_data_map_->at(target_name ) = target;
 
#ifdef __DEBUG_TENSOP_COMPUTER
  cout << "(post_summation) target->norm() = " << setprecision(13) <<  target->norm(); cout << endl << endl;
#endif
  return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Adds summand name factor to target name 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::add_tensors(string target_name, string summand_name, DataType factor ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::add_tensors " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   shared_ptr<Tensor_<DataType>> target_tens  = find_or_get_CTP_data( target_name );
   shared_ptr<Tensor_<DataType>> summand_tens = find_or_get_CTP_data( summand_name );
   Tensor_Calc_->add_tensors( target_tens, summand_tens, factor ); 

   return; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor T3 with elements T3_{ijkl..} = T1_{ijkl..}/T2_{ijkl..} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> TensOp_Computer::TensOp_Computer<DataType>::divide_tensors(string T1_name, string T2_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::divide_tensors " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
   return Tensor_Calc_->divide_tensors( find_or_get_CTP_data(T1_name),  find_or_get_CTP_data(T2_name)); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Builds a tensor with the relevant ranges and puts it into the map at the name
//Intended for intermediate and temporary tensors
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> TensOp_Computer::TensOp_Computer<DataType>::get_tensor( const std::vector<string>& id_ranges, bool set_elems, 
                                                                                      DataType init_val ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::get_tensor " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> new_tens; 
  if ( id_ranges.size() != 0 ) {
    new_tens = make_shared<Tensor_<DataType>>( Get_Bagel_IndexRanges( id_ranges ) );
  } else {
    new_tens = make_shared<Tensor_<DataType>>( vector<IndexRange>( 1, IndexRange(1,1,0,1) ) );
  }  
  new_tens->allocate();
  
  if ( set_elems ) {
    if ( init_val == (DataType)(0.0) ) {
      new_tens->zero(); 
    } else {
      Tensor_Arithmetic::Tensor_Arithmetic<DataType>::set_tensor_elems( new_tens, init_val );
    }  
  }
 
  return new_tens;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Builds a tensor with the relevant ranges and puts it into the map at the name
//Intended for intermediate and temporary tensors
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::build_tensor( string new_tens_name , const std::vector<string>& id_ranges, DataType init_val ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::build_tensor " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> new_tens = get_tensor( id_ranges, true /* = set_elems */, init_val );

  auto new_tens_loc = tensop_data_map_->find( new_tens_name); 
  if(  new_tens_loc != tensop_data_map_->end()){
    cout << "The tensor block " << new_tens_name << "is already in the map.... overwriting " << endl;
    new_tens_loc->second = new_tens;
  } else { 
    tensop_data_map_->emplace ( new_tens_name , new_tens ); 
  } 

  cout << "new_tens_name->norm() = "  << tensop_data_map_->at(new_tens_name) << endl;
 
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//As above, but does division in place, modifying the input;  T1_{ijkl..} = T1_{ijkl..}/T2_{ijkl..} 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::divide_tensors_in_place(string T1_name, string T2_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << " TensOp_Computer::TensOp_Computer::divide_tensors_in_place" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  Tensor_Calc_->divide_tensors_in_place( find_or_get_CTP_data(T1_name),  find_or_get_CTP_data(T2_name)); 

  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>> TensOp_Computer::TensOp_Computer<DataType>::get_block_Tensor(string tens_block_name){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::TensOp_Computer::get_block_Tensor : " << tens_block_name << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<Tensor_<DataType>> tens; 
  
  if(  tensop_data_map_->find(tens_block_name) != tensop_data_map_->end()){
    tens = tensop_data_map_->at(tens_block_name);

  } else {
    assert( false );
    cout <<"not in map ... " <<  tens_block_name << " must be formed from direct product tensor " << endl; 
    vector<string> sub_tens_names(0);
    for ( shared_ptr<CtrTensorPart_Base>& ctp : *(CTP_map_->at(tens_block_name)->CTP_vec()))
      sub_tens_names.push_back(ctp->name());
    
    tens = direct_product_tensors( sub_tens_names ); 

  }
  return tens;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO this can end up copying a lot; ideally should call the integrals just for these blocks
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
void TensOp_Computer::TensOp_Computer<DataType>::get_tensor_data_blocks(shared_ptr<set<shared_ptr<Range_Block_Info>>> required_blocks ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::TensOp_Computer::get_tensor_data_blocks " << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   for ( auto& block : *required_blocks ) { 
 
     string block_name = block->name();
     cout << "fetching block " << block_name ; cout.flush(); 
     
     if(  tensop_data_map_->find(block_name) != tensop_data_map_->end()){
        cout << " .. already in map" <<  endl;

     } else {

       string full_tens_name = block->op_info_->op_state_name_;

       shared_ptr<vector<string>> id_ranges = CTP_map_->at(block_name)->unc_id_ranges() ;
     
       if( tensop_data_map_->find(full_tens_name) == tensop_data_map_->end()){

         // TODO This will get the whole tensor, really, we should just get the blocks we want
         if ( full_tens_name[0] == 'H' || full_tens_name[0] == 'h' || full_tens_name[0] == 'f' ) {  
           build_mo_tensor( full_tens_name ); 
           get_sub_tensor( full_tens_name, block_name, *id_ranges );

         } else if ( full_tens_name[0] == 'X' || full_tens_name[0] == 'T' || full_tens_name[0] == 't'  ) {  
           build_tensor( block_name, *id_ranges, (DataType)(1.0) );
 //          shared_ptr<Tensor_<DataType>> test_block = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::get_test_tensor_row_major( tensop_data_map_->at(block_name)->indexrange() ); 
 //          tensop_data_map_->at(block_name) = test_block;
            
         } else if ( full_tens_name[0] == 'S') {  
          vector<string> test_block = { "c", "a", "v", "v" };
//          vector<string> test_block = { "a", "c", "v", "v" };
          // vector<string> test_block = { "c", "c", "v", "v" };
           if ( *id_ranges == test_block ){ 
             cout << "setting S QQQQQQQQQQQQQQ" << endl;
             build_tensor( block_name, *id_ranges, (DataType)(1.0) );
           } else { 
             build_tensor( block_name, *id_ranges, (DataType)(0.0) );
           }
         
         }
 
       } else {
       
         get_sub_tensor( full_tens_name, block_name, *id_ranges );
  
       }
     }
   } 
   return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Gets ranges and factors from the input which will be used in definition of terms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp_Computer::TensOp_Computer<DataType>::build_mo_tensor( string mo_tensor_name ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::TensOp_Computer::calculate_mo_integrals() : " << mo_tensor_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if ( mo_tensor_name[0] == 'H' ) {   
    auto H_loc = tensop_data_map_->find( "H_{00}" );

    if ( H_loc == tensop_data_map_->end() ) {
      if ( ! moint_computer_->got_fock_coeffs() ) {
        vector<string> free2 = { "free", "free" };
        moint_computer_->calculate_fock( free2, true, true );
        tensop_data_map_->emplace( "f_{00}" , moint_computer_->f1() );
      }

      vector<string> free4 = { "free", "free", "free", "free" };
      moint_computer_->calculate_v2( free4 ) ;
      tensop_data_map_->emplace( "H_{00}" , moint_computer_->v2() );
 
      if ( mo_tensor_name != "H_{00}" ) 
        tensop_data_map_->emplace( mo_tensor_name , moint_computer_->v2() );

    } else { 
      tensop_data_map_->emplace( mo_tensor_name , H_loc->second );

    }

  } else if ( mo_tensor_name[0] == 'h' ) {   
    auto f_loc = tensop_data_map_->find( "h_{00}" );

    if ( f_loc == tensop_data_map_->end() ) {
      vector<string> free2 = { "free", "free" };
      moint_computer_->calculate_h1( free2, false );
      tensop_data_map_->emplace( "h_{00}" , moint_computer_->h1() );
 
      if ( mo_tensor_name != "h_{00}" )
        tensop_data_map_->emplace( mo_tensor_name , moint_computer_->h1() );
    }

  } else if ( mo_tensor_name[0] == 'f' ) {   
    auto f_loc = tensop_data_map_->find( "f_{00}" );

    if ( f_loc == tensop_data_map_->end() ) {
      vector<string> free2 = { "free", "free" };
      moint_computer_->calculate_fock( free2, true, true );
      tensop_data_map_->emplace( "f_{00}" , moint_computer_->f1() );

      //TODO fix this properly for different states 
      if ( mo_tensor_name != "f_{00}" )
        tensop_data_map_->emplace( mo_tensor_name , moint_computer_->f1() );
    }
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_on_same_tensor( std::string T_in_name, std::string T_out_name, pair<int,int> ctr_todo) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
   cout << "TensOp_Computer::TensOp_Computer::contract_on_same_tensor "; cout.flush();
   cout << ": "  << T_in_name << " over (" << ctr_todo.first << ", " << ctr_todo.second << ") to get " << T_out_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   return Tensor_Calc_->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), ctr_todo  ); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_on_same_tensor( string T_in_name, shared_ptr<vector<int>> ctrs_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::TensOp_Computer::contract_on_same_tensor" << ": "  << T_in_name << " over "; print_vector( *ctrs_pos) ; cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   return Tensor_Calc_->contract_on_same_tensor( find_or_get_CTP_data(T_in_name), *ctrs_pos  ); 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType
TensOp_Computer::TensOp_Computer<DataType>::sum_elems( std::string Tens_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::dot_different_tensors" << endl;
cout << ": "  << Tens_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<Tensor_<DataType>> Tens = find_or_get_CTP_data(Tens_name);

  return Tens->norm();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
DataType
TensOp_Computer::TensOp_Computer<DataType>::dot_tensors( std::string T1_name, std::string T2_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::dot_different_tensors" << endl;
cout << ": "  << T1_name << " dot " << T2_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<Tensor_<DataType>> T1 = find_or_get_CTP_data(T1_name);
  shared_ptr<Tensor_<DataType>> T2 = find_or_get_CTP_data(T2_name);

  return T1->dot_product( T2 );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                                        pair<int,int> ctr_todo ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::contract_on_different_tensor" << endl;
cout << ": "  << T1_in_name << " and " << T2_in_name << " over T1[" << ctr_todo.first << "] and T2[" << ctr_todo.second << "] to get " << T_out_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<Tensor_<DataType>> Tens1_in = find_or_get_CTP_data(T1_in_name);
  shared_ptr<Tensor_<DataType>> Tens2_in = find_or_get_CTP_data(T2_in_name);
  return Tensor_Calc_->contract_different_tensors(Tens1_in, Tens2_in, ctr_todo );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::contract_different_tensors( std::string T1_in_name, std::string T2_in_name, std::string T_out_name,
                                                                        std::pair<std::vector<int>,std::vector<int>> ctrs_todo                  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::contract_on_different_tensor" << endl; 
cout << "T1_in_name = " << T1_in_name; cout.flush();
cout << "  T2_in_name = " << T2_in_name; cout.flush();
cout << "  T_out_name = " << T_out_name << endl;
cout << "contractions : [ " ; cout.flush(); 
for ( int qq = 0; qq != ctrs_todo.first.size(); qq++ )
  cout << "[ " << ctrs_todo.first[qq] << ", " << ctrs_todo.second[qq] << " ]" << endl;
cout << " ] " << endl;
assert( ctrs_todo.first.size() != ctrs_todo.second.size() );
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> Tens1_in = find_or_get_CTP_data(T1_in_name);
  shared_ptr<Tensor_<DataType>> Tens2_in = find_or_get_CTP_data(T2_in_name);

  return Tensor_Calc_->contract_different_tensors(Tens1_in, Tens2_in, ctrs_todo );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::direct_product_tensors( std::vector<std::string>& tensor_names ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::direct_product_tensors" << endl; 
print_vector(tensor_names, "tensor_names"); cout <<endl; 
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  assert(tensor_names.size() > 1 );

  string Tname_comp = tensor_names.front();
  shared_ptr<Tensor_<DataType>> Tens_prod = find_or_get_CTP_data( Tname_comp );
  shared_ptr<Tensor_<DataType>> Tens_intermediate;
  for ( vector<string>::iterator tn_it = tensor_names.begin()+1 ; tn_it != tensor_names.end(); tn_it++ ) {
    shared_ptr<Tensor_<DataType>> Tens_next = find_or_get_CTP_data(*tn_it);
    Tens_intermediate = Tensor_Arithmetic::Tensor_Arithmetic<DataType>::direct_tensor_product( Tens_prod, Tens_next ); 
    Tens_prod = std::move(Tens_intermediate);
    Tname_comp += *tn_it;
  }           
 
  return Tens_prod;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
shared_ptr<Tensor_<DataType>>
TensOp_Computer::TensOp_Computer<DataType>::reorder_block_Tensor(string tens_block_name, vector<int>& new_order){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER_VERBOSE
 cout << "TensOp_Computer::TensOp_Computer::reorder_block_Tensor "; cout.flush();
 cout << " : " << tens_block_name ; cout.flush();  WickUtils::print_vector( new_order, "    new_order"); cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  auto tensop_data_map_loc = tensop_data_map_->find( tens_block_name ); 
  shared_ptr<Tensor_<DataType>> T_part; 
  if( tensop_data_map_loc == tensop_data_map_->end() ){
    T_part =  get_block_Tensor( tens_block_name ); 
    tensop_data_map_->emplace( tens_block_name , T_part );
  } else { 
    T_part = tensop_data_map_loc->second; 
  }

  return Tensor_Calc_->reorder_block_Tensor( T_part , new_order );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
pair<int,int>
TensOp_Computer::TensOp_Computer<DataType>::relativize_ctr_positions( pair <int,int> ctr_todo, 
                                                                      shared_ptr<CtrTensorPart_Base> ctp1,
                                                                      shared_ptr<CtrTensorPart_Base> ctp2 ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
 cout << "TensOp_Computer::TensOp_Computer::relativize_ctr_positions" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////

   pair<int,int> rel_ctr;

   int T1_orig_size = ctp1->full_id_ranges()->size();
   int T2_orig_size = ctp2->full_id_ranges()->size();

   if (ctr_todo.first >= T1_orig_size ){
     rel_ctr = make_pair(ctr_todo.second, ctr_todo.first-T2_orig_size);
   } else {
     rel_ctr = make_pair(ctr_todo.first, ctr_todo.second-T1_orig_size);
   }

  for ( int ii = 0 ; ii != ctp1->unc_pos()->size() ; ii++ )
     if ( ctp1->unc_pos(ii) == rel_ctr.first){
       rel_ctr.first = ii;
       break;
     }

  for ( int ii = 0 ; ii != ctp2->unc_pos()->size() ; ii++ )
     if (ctp2->unc_pos(ii) == rel_ctr.second){ 
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
#ifdef __DEBUG_TENSOP_COMPUTER
 cout << "TensOp_Computer::TensOp_Computer::find_or_get_CTP_data  : " <<  CTP_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> tensor_block;
  auto data_loc =  tensop_data_map_->find(CTP_name); 
  if ( data_loc == tensop_data_map_->end() ){
    tensor_block = get_block_Tensor(CTP_name);   
    tensop_data_map_->emplace(CTP_name, tensor_block);   
  } else {
    tensor_block = data_loc->second;
  }
  
  return tensor_block;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
vector<IndexRange>
TensOp_Computer::TensOp_Computer<DataType>::Get_Bagel_IndexRanges( const vector<string>& ranges_str ){ 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::Get_Bagel_IndexRanges 1arg "; print_vector(ranges_str, "ranges_str" ) ; cout << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<IndexRange> ranges_Bagel(ranges_str.size());
  for ( int ii = 0 ; ii != ranges_str.size(); ii++)
    ranges_Bagel[ii] = *range_conversion_map_->at(ranges_str[ii]);

  return ranges_Bagel;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void TensOp_Computer::TensOp_Computer<DataType>::build_test_tensor( string test_tensor_name, vector<size_t> dimensions, vector<size_t> max_block_sizes ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_TENSOP_COMPUTER
cout << "TensOp_Computer::TensOp_Computer::build_test_tensor : " << test_tensor_name <<  endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto tt_loc = tensop_data_map_->find( test_tensor_name );
  if ( tt_loc != tensop_data_map_->end() )
    throw logic_error("already a tensor with named \"" + test_tensor_name + "\" in the map! Aborting!" ); 

  vector<IndexRange> index_ranges(dimensions.size()); 
  for ( int ii = 0 ; ii != dimensions.size(); ii++ ) 
    index_ranges[ii] = SMITH::IndexRange(dimensions[ii], max_block_sizes[ii] ); 

  shared_ptr<Tensor_<DataType>> test_tens = make_shared<Tensor_<DataType>>( index_ranges );
  test_tens->allocate();

  print_vector( max_block_sizes, "max_block_sizes"); cout.flush();   print_vector( dimensions, "   dimensions"); cout << endl;
  set_test_elems( test_tens, test_tensor_name  );

  tensop_data_map_->emplace( test_tensor_name, test_tens );

  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename DataType>
void
TensOp_Computer::TensOp_Computer<DataType>::get_sub_tensor( std::string full_tens_name, std::string block_name,
                                                            const vector<string>& range_names ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_TENSOP_COMPUTER
cout << "Tensop_Computer::get_sub_tensor" << endl;
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<DataType>> full_tens = find_or_get_CTP_data( full_tens_name );

  vector<IndexRange> block = Get_Bagel_IndexRanges( range_names ); 

  tensop_data_map_->emplace( block_name , Tensor_Arithmetic_Utils::get_sub_tensor( full_tens, block ) );
  
  return;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template < typename DataType >
void TensOp_Computer::TensOp_Computer<DataType>::print_compute_list(string output_name ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_TENSOP_COMPUTER
cout << "TensOp_Computer::print_contraction_list" << endl;  
#endif ///////////////////////////////////////////////////////////////////////////////////////////////////////

  auto compute_map_loc = ACompute_map_->find(output_name); 
  if ( compute_map_loc == ACompute_map_->end() )
    throw logic_error("compute list for " + output_name + " not in the map" ); 
  
  shared_ptr<vector<shared_ptr<CtrOp_base>>> compute_list = compute_map_loc->second; 

  cout << "=========================================================================================================" << endl;
  cout << output_name << endl;
  cout << "=========================================================================================================" << endl;
  for (shared_ptr<CtrOp_base> ctr_op : *compute_list ){
    if ( ctr_op->ctr_type()[0] == 'd' ){
      cout << "[" << ctr_op->T1name() << " , " << ctr_op->T2name() << " , (";
      cout << ctr_op->T1_ctr_abs_pos() << "," <<  ctr_op->T2_ctr_abs_pos() << ")" << " , " << ctr_op->Tout_name() << " ] " ; cout << ctr_op->ctr_type() << endl;
  
    } else if (ctr_op->ctr_type()[0] == 's' ){
      cout << "[" << ctr_op->T1name() << " , " << ctr_op->T1name() << " , (";
      cout << ctr_op->ctr_abs_pos().first << "," <<  ctr_op->ctr_abs_pos().second << ")" << " , " << ctr_op->Tout_name() << " ] " ;   cout << ctr_op->ctr_type()  << endl;

    } else { 
      cout << ctr_op->ctr_type() << " " << ctr_op->Tout_name() <<  endl;
    }
    
  }
  cout << "=========================================================================================================" << endl;
 
  return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class TensOp_Computer::TensOp_Computer<double>;
template class TensOp_Computer::TensOp_Computer<complex<double>>;
///////////////////////////////////////////////////////////////////////////////////////////////////////
