#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/smith/wicktool/equation_computer.h>
#include <src/util/f77.h>
using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::Tensor_Sorter;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Equation_Computer::Equation_Computer::Equation_Computer( std::shared_ptr<const SMITH_Info<double>> ref, std::shared_ptr<Equation<double>> eqn_info_in,
                                                         std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_in,
                                                         std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> CTP_data_map_in,
                                                         std::shared_ptr<std::map<std::string, std::shared_ptr<Tensor_<double>>>> Gamma_data_map_in ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  eqn_info =  eqn_info_in;
  nelea_ = ref->ciwfn()->det()->nelea();
  neleb_ = ref->ciwfn()->det()->neleb();
  ncore_ = ref->ciwfn()->ncore();
  norb_  = ref->ciwfn()->nact();
  nstate_ = ref->ciwfn()->nstates();
  cc_ = ref->ciwfn()->civectors();
  maxtile = ref->maxtile();
  
  GammaMap = eqn_info->GammaMap;
  CTP_map  = eqn_info->CTP_map;
  
  CTP_data_map  = CTP_data_map_in;
  Gamma_data_map = Gamma_data_map_in;
  
  CIvec_data_map = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();
  Sigma_data_map = make_shared<std::map<std::string, std::shared_ptr<Tensor_<double>>>>();
  determinants_map = make_shared<std::map<std::string, std::shared_ptr<const Determinants>>>();

  cimaxblock = 2500; //figure out what is best, this chosen so 4 orbital indexes to hit maxtile of 10000.

  get_civector_indexranges(1);
  range_conversion_map = range_conversion_map_in;

  tester();
    
}  

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string Equation_Computer::Equation_Computer::get_det_name(shared_ptr<const Determinants> Detspace ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 return "[" +to_string(Detspace->norb()) + ",{"+to_string(Detspace->nelea())+"a,"+to_string(Detspace->neleb())+"b}]";
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//should be extended to deal with spin sectors
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::get_civector_indexranges(int nstates) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
  for ( int ii = 0 ; ii != nstates; ii++ ) 
    range_conversion_map->emplace( get_civec_name( ii , cc_->data(ii)->det()->norb(), cc_->data(ii)->det()->nelea(), cc_->data(ii)->det()->neleb()),
                                   make_shared<IndexRange>(cc_->data(ii)->lena()*cc_->data(ii)->lenb(), cimaxblock ));  


  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::tester(){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  shared_ptr<Tensor_<double>> civec1 =  convert_civec_to_tensor( cc_->data(0), 0 );
  shared_ptr<Tensor_<double>> civec2 =  convert_civec_to_tensor( cc_->data(0), 0 );

  double normval = civec1->dot_product(civec2); 
  cout << " civec1->dot_product(civec2) = " << normval << endl;
  cout << " civec1->rms()               = " << civec1->rms()  << endl;
  cout << " civec1->norm()              = " << civec1->norm() << endl;
  
  assert(!(abs(normval -1.00) > 0.000000000001) ); 
  
  return;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::Calculate_CTP(std::string A_contrib ){
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << " Equation_Computer::Equation_Computer::Calculate_CTP(std::string A_contrib_name )" << endl;
  if (eqn_info->ACompute_map->at(A_contrib)->size() == 0 )
    cout << "THIS COMPUTE LIST IS EMPTY" << endl;

  for (shared_ptr<CtrOp_base> ctr_op : *(eqn_info->ACompute_map->at(A_contrib))){
 
    // check if this is an uncontracted multitensor (0,0) && check if the data is in the map
    if( CTP_data_map->find(ctr_op->Tout_name()) == CTP_data_map->end() ) {
       shared_ptr<Tensor_<double>>  New_Tdata; 
 
      if ( ctr_op->ctr_type()[0] == 'g'){  cout << " : no contraction, fetch this tensor part" << endl; 
        New_Tdata =  get_block_Tensor(ctr_op->Tout_name());
        CTP_data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
  
      } else if ( ctr_op->ctr_type()[0] == 'd' ){ cout << " : contract different tensors" << endl; 
        New_Tdata = contract_different_tensors( make_pair(ctr_op->T1_ctr_rel_pos(), ctr_op->T2_ctr_rel_pos()), ctr_op->T1name(), ctr_op->T2name(), ctr_op->Tout_name());
        CTP_data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
      
      } else if ( ctr_op->ctr_type()[0] == 's' ) { cout << " : contract on same tensor" <<  endl; 
        New_Tdata = contract_on_same_tensor( ctr_op->ctr_rel_pos(), ctr_op->T1name(), ctr_op->Tout_name()); 
        CTP_data_map->emplace(ctr_op->Tout_name(), New_Tdata); 
      } else { 
        throw std::runtime_error(" unknown contraction type : " + ctr_op->ctr_type() ) ;
      }
    }
    cout << "A_contrib.first = " << A_contrib << endl;
    cout << "CTPout_name =  " << ctr_op->Tout_name() << endl;
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_on_same_tensor( pair<int,int> ctr_todo, std::string Tname, std::string Tout_name) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::contract_on_same_tensor" << endl;

   shared_ptr<CtrTensorPart<double>> CTP_old = CTP_map->at(Tname);
   shared_ptr<Tensor_<double>> CTP_data_old = find_or_get_CTP_data(Tname); cout <<  endl;

   // get original uncontracted ranges and positions of Ctrs relative to the current tensor
   int ctr1 = ctr_todo.first;
   int ctr2 = ctr_todo.second;
   vector<IndexRange> unc_ranges_old = CTP_data_old->indexrange(); 
   vector<IndexRange> unc_ranges_new(unc_ranges_old.size()-2);  
   vector<int> unc_pos_new(unc_ranges_old.size()-2);  
         
   vector<IndexRange>::iterator urn_iter = unc_ranges_new.begin();
   vector<int>::iterator upn_iter = unc_pos_new.begin();
   for ( int ii = 0 ; ii != CTP_old->unc_pos->size() ; ii++ )
     if ( (ii != ctr_todo.first) && (ii != ctr_todo.second) ){ 
       *urn_iter++ = unc_ranges_old[ii];
       *upn_iter++ = ii;
     }

   shared_ptr<Tensor_<double>> CTP_data_out = make_shared<Tensor_<double>>(unc_ranges_new);
   CTP_data_out->allocate();
   int num_ctr_blocks = unc_ranges_old[ctr_todo.first].range().size();

   //loops over index blocks where ctr1 = ctr2 
   for (int ii = 0 ; ii != num_ctr_blocks ; ii++){ 
     shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> maxs = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     for ( int jj = 0 ; jj != unc_ranges_old.size() ; jj++ ) 
        maxs->at(jj) = unc_ranges_old[jj].range().size()-1;
 
     mins->at(ctr_todo.first)  = ii;
     mins->at(ctr_todo.second) = ii; 
     maxs->at(ctr_todo.first)  = ii;
     maxs->at(ctr_todo.second) = ii;
     block_pos->at(ctr_todo.first)  = ii;
     block_pos->at(ctr_todo.second) = ii;
     
     const int ione = 1; 
     const double done = 1.0; 
     do {
       
       shared_ptr<vector<pair<size_t, size_t>>> bob = get_block_start( make_shared<vector<IndexRange>>(unc_ranges_old), block_pos ) ;
       vector<Index> CTP_id_blocks_old = *(get_rng_blocks(block_pos, unc_ranges_old)); 
       vector<Index> CTP_id_blocks_new(CTP_id_blocks_old.size()-2);
       for (int kk = 0 ; kk != unc_pos_new.size(); kk++)       
         CTP_id_blocks_new[kk] = CTP_id_blocks_old[unc_pos_new[kk]];
     
       vector<int> range_sizes = get_sizes(CTP_id_blocks_old);
       int total_size = accumulate( range_sizes.begin(), range_sizes.end(), 1, std::multiplies<int>() );

       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(range_sizes);
       int inner_stride = Tens_strides->at(ctr1) < Tens_strides->at(ctr2) ? Tens_strides->at(ctr1) : Tens_strides->at(ctr2);

       shared_ptr<vector<int>> maxs2 = make_shared<vector<int>>(range_sizes.size(),0);
       shared_ptr<vector<int>> mins2 = make_shared<vector<int>>(maxs->size(), 0); 
       shared_ptr<vector<int>> fvec2 = make_shared<vector<int>>(*mins); 

       //within index block, loop through data copying chunks where ctr1 = ctr2 
       //Has odd striding; probably very inefficient when ctr1 or ctr2 are the leading index.
       for (int jj = ctr1+1 ; jj != maxs2->size(); jj++) 
         maxs2->at(jj) = range_sizes[jj]-1; // note, copying inner block as vector, so want max = min = 0 for those indexes;
       
       int ctr1_rlen = range_sizes[ctr1];          
       int tmp = total_size/(ctr1_rlen*ctr1_rlen); 
       unique_ptr<double[]>      CTP_data_block_old = CTP_data_old->get_block(CTP_id_blocks_old);
       std::unique_ptr<double[]> CTP_data_block_new(new double[tmp]);
       std::fill_n(CTP_data_block_new.get(), tmp, 0.0);
       shared_ptr<vector<int>>   CTens_strides   =  get_CTens_strides( range_sizes, ctr1, ctr2 );
       
       for ( int jj = 0 ; jj != ctr1_rlen ; jj++ ) {

         fill(fvec2->begin(), fvec2->end(), 0);
         maxs2->at(ctr1) = jj; 
         maxs2->at(ctr2) = jj; 
         mins2->at(ctr1) = jj;  
         mins2->at(ctr2) = jj; 
         fvec2->at(ctr1) = jj; 
         fvec2->at(ctr2) = jj; 

         do { 
           int pos   = inner_product( fvec2->begin(), fvec2->end(), CTens_strides->begin(), 0); 
           int inpos = inner_product( fvec2->begin(), fvec2->end(), Tens_strides->begin(),  0); 

           blas::ax_plus_y_n(1.0, CTP_data_block_old.get()+inpos, inner_stride, CTP_data_block_new.get()+pos);

         } while (fvec_cycle_skipper_f2b(fvec2, maxs2, mins2));
       }
 
       CTP_data_out->add_block( CTP_data_block_new, CTP_id_blocks_new );
     } while (fvec_cycle_skipper(block_pos, maxs, mins ));
  } 

  CTP_map->at(Tout_name)->got_data = true; 
  cout << "CTP_data_out->rms()  = " + Tout_name + "->rms() =" << CTP_data_out->rms() << endl;

  return CTP_data_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Will contract over all the indexes given in the conctracted indexes list. 
//Note these are indexes relative to the input tensor, not the original Tensor from which the CtrTensorPart was obtained.
//In principal this could replace the contract_same_tensor routine, as there's a lot of duplicated code
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_on_same_tensor( vector<int>& contracted_index_positions, std::shared_ptr<Tensor_<double>> Tens_in ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_trace_tensor" << endl;

   if ( contracted_index_positions.size()-  Tens_in->indexrange().size() )  {
     cout << "contracted_index_positions         = [ "; cout.flush();  for ( int pos : contracted_index_positions ) {cout << pos << " " ; } cout << " ] " << endl;
     cout << "contracted_index_positions.size()  = "<<  contracted_index_positions.size() << endl;
     cout << "Tens_in->indexrange().size()       = " << Tens_in->indexrange().size()  <<  endl;
     cout << "trying to contract more indexes than there are in the tensor!" << endl;
     assert(false);
   }

   int num_ctrs = contracted_index_positions.size();
   // get original uncontracted ranges and positions of Ctrs relative to the current tensor
   vector<IndexRange> unc_ranges_old = Tens_in->indexrange(); 
   vector<IndexRange> unc_ranges_new(unc_ranges_old.size() - num_ctrs);  
   vector<int> unc_pos_new(unc_ranges_old.size() - num_ctrs);  

   vector<bool> unc_get(Tens_in->indexrange().size(), true);
   for( int ii = 0; ii !=contracted_index_positions.size(); ii++ )
      unc_get[contracted_index_positions[ii]] = false;

   for ( int ii = 0 ; ii != unc_get.size(); ii++  )
     if (unc_get[ii]) {
       unc_ranges_new[ii] = unc_ranges_old[ii];
       unc_pos_new[ii] = ii;
     }

   shared_ptr<Tensor_<double>> Tens_out = make_shared<Tensor_<double>>(unc_ranges_new);
   Tens_out->allocate();
   Tens_out->zero();

   // loops over index blocks with same ctr
   int num_ctr_blocks = unc_ranges_old[contracted_index_positions[0]].range().size();
   for (int ii = 0 ; ii != num_ctr_blocks ; ii++){ 
     shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> mins = make_shared<vector<int>>(unc_ranges_old.size(),0);  
     shared_ptr<vector<int>> maxs = make_shared<vector<int>>(unc_ranges_old.size());  
     for ( int jj = 0 ; jj != unc_ranges_old.size() ; jj++ ) 
        maxs->at(jj) = unc_ranges_old[jj].range().size()-1;
 
     for ( int ctr_idx_pos : contracted_index_positions ) {
        mins->at(ctr_idx_pos)       = ii;
        maxs->at(ctr_idx_pos)       = ii;
        block_pos->at(ctr_idx_pos)  = ii;
     }
     
     do {
       
       vector<Index> Tens_id_blocks_old = *(get_rng_blocks(block_pos, unc_ranges_old)); 
       vector<Index> Tens_id_blocks_new(Tens_id_blocks_old.size()- num_ctrs );
       for (int kk = 0 ; kk != unc_pos_new.size(); kk++)       
         Tens_id_blocks_new[kk] = Tens_id_blocks_old[unc_pos_new[kk]];
    
 
       vector<int> range_sizes = get_sizes(Tens_id_blocks_old);
       shared_ptr<vector<int>> Tens_strides = get_Tens_strides(range_sizes);
   
       // To save time,  the contracted parts are copied in blocks i.e. (i, j, k, k, : , : ) 
       // Inner stride is the size of the data block to the right of the contracted indexes
       int inner_stride = Tens_strides->at(0); 
       int innermost_ctr = contracted_index_positions.back() ;
       for ( int qq = 0 ; qq!=contracted_index_positions.size() ; qq++ ){
         if ( (inner_stride) > Tens_strides->at(contracted_index_positions[qq])){
           inner_stride = Tens_strides->at(contracted_index_positions[qq]);
           innermost_ctr = contracted_index_positions[qq];
         }
       }

       cout << "inner_stride  = " << inner_stride << endl;
       cout << "innermost_ctr = " << innermost_ctr << endl;

       shared_ptr<vector<int>> mins2 = make_shared<vector<int>>(maxs->size(), 0); 
       shared_ptr<vector<int>> fvec2 = make_shared<vector<int>>(*mins); 
       shared_ptr<vector<int>> maxs2 = make_shared<vector<int>>(range_sizes.size(),0);
       // note: copying inner block as vector, so want max = min = 0 for those indexes  
       // looks backwards, but is countered by use of f2b fvec routine
       for (int jj = innermost_ctr ; jj != maxs2->size(); jj++) 
         maxs2->at(jj) = range_sizes[jj]-1; 

       shared_ptr<vector<int>> CTens_strides = get_CTens_strides( range_sizes, contracted_index_positions );

       int ctr_block_len = range_sizes[contracted_index_positions.front()];          
       int out_block_size = 1;
       for ( int pos : contracted_index_positions) 
          out_block_size *= range_sizes[pos];

       unique_ptr<double[]>      Tens_data_block_old = Tens_in->get_block(Tens_id_blocks_old);
       std::unique_ptr<double[]> Tens_data_block_new(new double[out_block_size]);
       std::fill_n(Tens_data_block_new.get(), out_block_size, 0.0);
       
       for ( int jj = 0 ; jj != ctr_block_len ; jj++ ) {

         fill(fvec2->begin(), fvec2->end(), 0);
         for ( int ctr_idx_pos : contracted_index_positions ) {
           maxs2->at(ctr_idx_pos) = jj; 
           mins2->at(ctr_idx_pos) = jj;  
           fvec2->at(ctr_idx_pos) = jj; 
         }

         //within index block, loop through data copying chunks where contracted_ctrs are equal
         do { 
           int pos   = inner_product( fvec2->begin(), fvec2->end(), CTens_strides->begin(), 0); 
           int inpos = inner_product( fvec2->begin(), fvec2->end(), Tens_strides->begin(),  0); 
           cout << "pos   = "<< pos <<   endl;
           cout << "inpos = "<< inpos <<   endl;
            
           blas::ax_plus_y_n(1.0, Tens_data_block_old.get()+inpos, inner_stride, Tens_data_block_new.get()+pos);
  
      
         } while (fvec_cycle_skipper_f2b(fvec2, maxs2, mins2));
       }
 
       Tens_out->add_block( Tens_data_block_new, Tens_id_blocks_new );
     } while (fvec_cycle_skipper(block_pos, maxs, mins ));
  } 

  return Tens_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Contracts tensors T1 and T2 over two specified indexes
//T1_org_rg T2_org_rg are the original index ranges for the tensors (not necessarily normal ordered).
//T2_new_rg T1_new_rg are the new ranges, with the contracted index at the end, and the rest in normal ordering.
//T1_new_order and T2_new_order are the new order of indexes, and are used for rearranging the tensor data.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>>
Equation_Computer::Equation_Computer::contract_different_tensors( pair<int,int> ctr_todo, std::string T1name, std::string T2name,
                                                                  std::string Tout_name){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::contract_on_different_tensor" <<endl; 

  shared_ptr<CtrTensorPart<double>> CTP1 = CTP_map->at(T1name); 
  shared_ptr<CtrTensorPart<double>> CTP2 = CTP_map->at(T2name); 

  pair<int,int> ctr_todo_rel = ctr_todo;
  shared_ptr<vector<int>> T1_org_order= make_shared<vector<int>>(CTP1->unc_pos->size());
  for (int ii =0 ; ii != T1_org_order->size(); ii++) { T1_org_order->at(ii) = ii ;}

  shared_ptr<vector<int>> T2_org_order= make_shared<vector<int>>(CTP2->unc_pos->size());
  for (int ii =0 ; ii != T2_org_order->size(); ii++) { T2_org_order->at(ii) = ii ;}

  shared_ptr<Tensor_<double>> CTP1_data = find_or_get_CTP_data(T1name);
  shared_ptr<vector<int>> T1_new_order  = put_ctr_at_back( T1_org_order , ctr_todo_rel.first);
  shared_ptr<vector<shared_ptr<const IndexRange>>> T1_org_rngs = Get_Bagel_const_IndexRanges(CTP1->full_id_ranges, CTP1->unc_pos) ;
  shared_ptr<vector<shared_ptr<const IndexRange>>> T1_new_rngs = reorder_vector(T1_new_order, T1_org_rngs);
  shared_ptr<vector<int>> maxs1 = get_num_index_blocks_vec(T1_new_rngs) ;

  shared_ptr<Tensor_<double>> CTP2_data = find_or_get_CTP_data(T2name);
  shared_ptr<vector<int>> T2_new_order  = put_ctr_at_front( T2_org_order , ctr_todo_rel.second);
  shared_ptr<vector<shared_ptr<const IndexRange>>> T2_org_rngs = Get_Bagel_const_IndexRanges(CTP2->full_id_ranges, CTP2->unc_pos) ;
  shared_ptr<vector<shared_ptr<const IndexRange>>> T2_new_rngs = reorder_vector(T2_new_order, T2_org_rngs);
  shared_ptr<vector<int>> maxs2 = get_num_index_blocks_vec(T2_new_rngs) ;

  auto Tout_unc_rngs = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T1_new_rngs->begin(), T1_new_rngs->end()-1);
  Tout_unc_rngs->insert(Tout_unc_rngs->end(), T2_new_rngs->begin()+1, T2_new_rngs->end());

  shared_ptr<Tensor_<double>> T_out;  

  {
  auto Tout_unc_rngs_raw = make_shared<vector<IndexRange>>(0); 
  for (auto id_rng : *Tout_unc_rngs)
    Tout_unc_rngs_raw->push_back(*id_rng);
  T_out = make_shared<Tensor_<double>>(*Tout_unc_rngs_raw); 
  }
  T_out->allocate();
  T_out->zero();

  //loops over all index blocks of T1 and  T2; final index of T1 is same as first index of T2 due to contraction
  auto T1_rng_block_pos = make_shared<vector<int>>(T1_new_order->size(),0);

  do { 
    cout << T1name << " block pos =  [ " ;    for (int block_num : *T1_rng_block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
    shared_ptr<vector<Index>> T1_new_rng_blocks = get_rng_blocks( T1_rng_block_pos, T1_new_rngs); 
    shared_ptr<vector<Index>> T1_org_rng_blocks = inverse_reorder_vector(T1_new_order, T1_new_rng_blocks); 
    
    size_t ctr_block_size = T1_new_rng_blocks->back().size(); 
    size_t T1_unc_block_size = get_block_size( T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1); 
    size_t T1_block_size =  get_block_size(T1_org_rng_blocks->begin(), T1_org_rng_blocks->end()); 

    std::unique_ptr<double[]> T1_data_new;
    {
      std::unique_ptr<double[]> T1_data_org = CTP1_data->get_block(*T1_org_rng_blocks);
      T1_data_new = reorder_tensor_data(T1_data_org.get(), T1_new_order, T1_org_rng_blocks);
    }
    
    shared_ptr<vector<int>> T2_rng_block_pos = make_shared<vector<int>>(T2_new_order->size(), 0);
    T2_rng_block_pos->front() = T1_rng_block_pos->back();
 
    do { 

      cout << T2name << " block pos =  [ " ;    for (int block_num : *T2_rng_block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
      shared_ptr<vector<Index>> T2_new_rng_blocks = get_rng_blocks(T2_rng_block_pos, T2_new_rngs); 
      shared_ptr<vector<Index>> T2_org_rng_blocks = inverse_reorder_vector(T2_new_order, T2_new_rng_blocks); 
      size_t T2_unc_block_size = get_block_size(T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      std::unique_ptr<double[]> T2_data_new;
      {
        std::unique_ptr<double[]> T2_data_org = CTP2_data->get_block(*T2_org_rng_blocks);
        T2_data_new = reorder_tensor_data(T2_data_org.get(), T2_new_order, T2_org_rng_blocks); 
      }
       
      std::unique_ptr<double[]> T_out_data(new double[T1_unc_block_size*T2_unc_block_size]);
      std::fill_n(T_out_data.get(), T1_unc_block_size*T2_unc_block_size, 0.0);
      //should not use transpose; instead build T2_new_order backwards... 
      dgemm_("N", "N", T1_unc_block_size, T2_unc_block_size, ctr_block_size, 1.0, T1_data_new, T1_unc_block_size,
              T2_data_new, ctr_block_size, 0.0, T_out_data, T1_unc_block_size );

      vector<Index> T_out_rng_block(T1_new_rng_blocks->begin(), T1_new_rng_blocks->end()-1);
      T_out_rng_block.insert(T_out_rng_block.end(), T2_new_rng_blocks->begin()+1, T2_new_rng_blocks->end());
      T_out->add_block( T_out_data, T_out_rng_block );

      //remove last index; contracted index is cycled in T1 loop
    } while(fvec_cycle_skipper(T2_rng_block_pos, maxs2, 0 )) ;

  } while (fvec_cycle_test(T1_rng_block_pos, maxs1 ));

  cout << "CTP1_data->rms()  = " + T1name + "->rms() =" << CTP1_data->rms() << endl;
  cout << "CTP2_data->rms()  = " + T2name + "->rms() =" << CTP2_data->rms() << endl;
  cout << "CTP_data_out->rms()  = " + Tout_name + "->rms() =" << T_out->rms() << endl;
  CTP_map->at(Tout_name)->got_data = true; 
  return T_out;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Sets all elements of input tensor Tens to the value specified by elem_val 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Equation_Computer::Equation_Computer::set_tensor_elems(shared_ptr<Tensor_<double>> Tens, double elem_val  ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::set_tensor_elems  " << endl;
  
   vector<IndexRange> id_ranges =  Tens->indexrange();

   shared_ptr<vector<int>> range_lengths = make_shared<vector<int>>(0); 
   for (IndexRange idrng : id_ranges )
     range_lengths->push_back(idrng.range().size()-1); 

//   cout << " range_lengths = [ " ; cout.flush(); for ( int pos : *range_lengths ) { cout << pos << " "  ; cout.flush(); } cout << "] " << endl;

   shared_ptr<vector<int>> block_pos = make_shared<vector<int>>(range_lengths->size(),0);  
   shared_ptr<vector<int>> mins = make_shared<vector<int>>(range_lengths->size(),0);  
   do {
     // cout << " block_pos = [ " ; cout.flush(); for ( int pos : *block_pos ) { cout << pos << " "  ; cout.flush(); } cout << "] " << endl;
     vector<Index> id_blocks(id_ranges.size());
     for( int ii = 0 ;  ii != id_blocks.size(); ii++)
       id_blocks[ii] =  id_ranges[ii].range(block_pos->at(ii));
 
     if ( Tens->exists(id_blocks) ) { 
       unique_ptr<double[]> block_data = Tens->get_block(id_blocks);
       std::fill_n(block_data.get(), Tens->get_size(id_blocks), elem_val);
       Tens->put_block(block_data, id_blocks);
     }
   } while (fvec_cycle_skipper(block_pos, range_lengths, mins ));

   return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_block_Tensor(string Tname){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_block_Tensor : " << Tname << endl;
  
   if(  CTP_data_map->find(Tname) != CTP_data_map->end())
     return CTP_data_map->at(Tname);

   shared_ptr<vector<string>> unc_ranges = CTP_map->at(Tname)->unc_id_ranges;  

   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   shared_ptr<vector<int>> range_lengths = make_shared<vector<int>>(0); 
   for (auto idrng : *Bagel_id_ranges )
     range_lengths->push_back(idrng.range().size()-1); 

   shared_ptr<Tensor_<double>> fulltens = CTP_data_map->at(Tname.substr(0,1));
   shared_ptr<Tensor_<double>> block_tensor = make_shared<Tensor_<double>>(*Bagel_id_ranges);
   block_tensor->allocate();
   block_tensor->zero();

   auto block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   auto mins = make_shared<vector<int>>(unc_ranges->size(),0);  
   do {
     cout << Tname << " block pos =  [ " ;    for (int block_num : *block_pos )  { cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
     
     vector<Index> T_id_blocks(Bagel_id_ranges->size());
     for( int ii = 0 ;  ii != T_id_blocks.size(); ii++)
       T_id_blocks[ii] =  Bagel_id_ranges->at(ii).range(block_pos->at(ii));
  
     unique_ptr<double[]> T_block_data = fulltens->get_block(T_id_blocks);
     block_tensor->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));

   return block_tensor;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a block of a tensor, defined as a new tensor, is copying needlessly, so find another way. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::reorder_block_Tensor(string Tname, shared_ptr<vector<int>> new_order){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::reorder_block_Tensor "; cout.flush();
   cout << " : " << Tname ; cout.flush();
   cout <<  " New_order = [ "; cout.flush();  for (int pos : *new_order ) { cout << pos << " " ; cout.flush(); } cout << "] " << endl;
  
   shared_ptr<Tensor_<double>> orig_tens ;
   { 
   auto CTP_data_map_loc = CTP_data_map->find(Tname); 
   if(  CTP_data_map_loc == CTP_data_map->end()){
     throw std::runtime_error(" don't have tensor data for " +Tname+ " yet.... Equation_Computer::reorder_block_Tensor " ) ;
   } else { 
     orig_tens =  CTP_data_map_loc->second; 
   } 
   }

   shared_ptr<vector<string>> unc_ranges = CTP_map->at(Tname)->unc_id_ranges;  
   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);
   shared_ptr<vector<int>> range_lengths = make_shared<vector<int>>(0); 
   for (auto idrng : *Bagel_id_ranges )
     range_lengths->push_back(idrng.range().size()-1); 
   
   shared_ptr<vector<IndexRange>> reordered_ranges = reorder_vector(new_order, Bagel_id_ranges ) ;
   shared_ptr<Tensor_<double>> reordered_block_tensor = make_shared<Tensor_<double>>(*reordered_ranges);
   reordered_block_tensor->allocate();
   reordered_block_tensor->zero();

   auto block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   auto mins = make_shared<vector<int>>(unc_ranges->size(),0);  
   do {
     cout << Tname << " block pos =  [ " ; for(int block_num : *block_pos ){ cout << block_num << " " ; cout.flush(); } cout << " ] " << endl;
     
     shared_ptr<vector<Index>>  orig_id_blocks = make_shared<vector<Index>>(Bagel_id_ranges->size());
     for( int ii = 0 ;  ii != orig_id_blocks->size(); ii++)
       orig_id_blocks->at(ii) =  Bagel_id_ranges->at(ii).range(block_pos->at(ii));
  
     unique_ptr<double[]> orig_data_block = orig_tens->get_block(*orig_id_blocks);
     shared_ptr<vector<Index>> reordered_id_blocks = reorder_vector(new_order, orig_id_blocks ) ;

     unique_ptr<double[]> reordered_data_block = reorder_tensor_data( orig_data_block.get(), new_order, orig_id_blocks ) ;
      
     reordered_block_tensor->put_block(reordered_data_block, *reordered_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));

   return reordered_block_tensor;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Returns a tensor with ranges specified by unc_ranges, where all values are equal to XX  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<Tensor_<double>> Equation_Computer::Equation_Computer::get_uniform_Tensor(shared_ptr<vector<string>> unc_ranges, double XX ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   cout << "Equation_Computer::Equation_Computer::get_uniform_Tensor" << endl;

   shared_ptr<vector<IndexRange>> Bagel_id_ranges = Get_Bagel_IndexRanges(unc_ranges);

   auto  range_lengths  = make_shared<vector<int>>(0); 
   for (auto idrng : *Bagel_id_ranges )
      range_lengths->push_back(idrng.range().size()-1); 

   shared_ptr<Tensor_<double>> block_tensor = make_shared<Tensor_<double>>(*Bagel_id_ranges);
   block_tensor->allocate();

   auto block_pos = make_shared<vector<int>>(unc_ranges->size(),0);  
   auto mins = make_shared<vector<int>>(unc_ranges->size(),0);  
   do {

     vector<Index> T_id_blocks(Bagel_id_ranges->size());
     for( int ii = 0 ;  ii != T_id_blocks.size(); ii++)
       T_id_blocks[ii] =  Bagel_id_ranges->at(ii).range(block_pos->at(ii));
     
     int out_size = 1;
     for ( Index id : T_id_blocks)
        out_size*= id.size();

     unique_ptr<double[]> T_block_data( new double[out_size] );
   
     double* dit = T_block_data.get() ;
     for ( int qq =0 ; qq != out_size; qq++ ) 
        T_block_data[qq] = XX; 

     block_tensor->put_block(T_block_data, T_id_blocks);

   } while (fvec_cycle(block_pos, range_lengths, mins ));
   return block_tensor;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
unique_ptr<double[]>
Equation_Computer::Equation_Computer::get_block_of_data( double* data_ptr ,
                                                         shared_ptr<vector<IndexRange>> id_ranges, 
                                                         shared_ptr<vector<int>> block_pos) {
////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Equation_Computer::get_block_of_data" << endl;

 // merge this to one loop, but keep until debugged, as this is likely to go wrong...
 // getting id position of id_block ; list of block sizes looks like [n,n,n,n,... , n-1, n-1 n-1] 
  vector<size_t> id_pos(block_pos->size());
  cout <<endl << "id_pos = ";
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){

    const size_t range_size         = id_ranges->at(ii).size();
    const size_t biggest_block_size = id_ranges->at(ii).range(0).size();
    const size_t num_blocks         = id_ranges->at(ii).range().size();
    const size_t remainder          = num_blocks * biggest_block_size - range_size;

    if (block_pos->at(ii) <= remainder  ){
       id_pos[ii] = num_blocks*block_pos->at(ii);//  + id_ranges->at(ii).range(block_pos->at(ii)).offset();

    } else if ( block_pos->at(ii) > remainder ) {
       id_pos[ii] = num_blocks*(range_size - remainder)+(num_blocks-1)*(remainder - block_pos->at(ii));// + id_ranges->at(ii).range(block_pos->at(ii)).offset(); 
    }; 
    cout << id_pos[ii] << " " ;
  }

  cout << endl << "range_sizes = " ;
  // getting size of ranges (seems to be correctly offset for node)
  vector<size_t> range_sizes(block_pos->size());
  for (int ii = 0 ; ii != id_ranges->size() ; ii++){
    range_sizes[ii]  = id_ranges->at(ii).size();
    cout << range_sizes[ii] << " " ;
  }

  size_t data_block_size = 1;
  size_t data_block_pos  = 1;
  for (int ii = 0 ; ii != id_ranges->size()-1 ; ii++){
    data_block_pos  *= id_pos[ii]*range_sizes[ii];
    data_block_size *= id_ranges->at(ii).range(block_pos->at(ii)).size();
  }

  data_block_size *= id_ranges->back().range(block_pos->back()).size();
 
  cout << "data_block_size = " << data_block_size << endl;
  cout << "data_block_pos = "  << data_block_pos << endl;
 
  unique_ptr<double[]> data_block(new double[data_block_size])  ;

  copy_n(data_block.get(), data_block_size, data_ptr+data_block_pos);

  return data_block; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<int,int>
Equation_Computer::Equation_Computer::relativize_ctr_positions(pair <int,int> ctr_todo, 
                                                               shared_ptr<CtrTensorPart<double>>  CTP1,
                                                               shared_ptr<CtrTensorPart<double>>  CTP2){
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
  auto CTP_data_loc =  CTP_data_map->find(CTP_name); 
  if ( CTP_data_loc == CTP_data_map->end() ){
    CTP_data = get_block_Tensor(CTP_name);   
  } else {
    CTP_data = CTP_data_loc->second;
  }
  
  return CTP_data;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
 Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 1arg" << endl;

  auto ranges_Bagel = make_shared<vector<shared_ptr<const IndexRange>>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(make_shared<const IndexRange>(*range_conversion_map->at(rng)));

  return ranges_Bagel;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<shared_ptr<const IndexRange>>>
Equation_Computer::Equation_Computer::Get_Bagel_const_IndexRanges(shared_ptr<vector<string>> ranges_str, shared_ptr<vector<int>> unc_pos){ 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_const_IndexRanges 2arg" << endl;

  vector<shared_ptr<const IndexRange>>  ranges_Bagel(unc_pos->size());
  for ( int ii = 0 ; ii != unc_pos->size() ; ii++) 
    ranges_Bagel[ii]=(make_shared<const IndexRange>(*range_conversion_map->at(ranges_str->at(unc_pos->at(ii)))));

  return make_shared<vector<shared_ptr<const IndexRange>>>(ranges_Bagel);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<IndexRange>> Equation_Computer::Equation_Computer::Get_Bagel_IndexRanges(shared_ptr<vector<string>> ranges_str){ 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "Get_Bagel_IndexRanges 1arg" << endl;

  auto ranges_Bagel = make_shared<vector<IndexRange>>(0);
  for ( auto rng : *ranges_str) 
    ranges_Bagel->push_back(*range_conversion_map->at(rng));

  return ranges_Bagel;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::inverse_reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "inverse_reorder_vector" << endl;

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
    newvec->at(neworder->at(ii)) =  origvec->at(ii);

  return newvec;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class vtype>
shared_ptr<vector<vtype>> Equation_Computer::Equation_Computer::reorder_vector(shared_ptr<vector<int>> neworder , shared_ptr<vector<vtype>> origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
cout << "reorder_vector" << endl;

  auto newvec = make_shared<vector<vtype>>(origvec->size());
  for( int ii = 0; ii != origvec->size(); ii++ )
     newvec->at(ii) = origvec->at(neworder->at(ii));
  return newvec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class DataType>
unique_ptr<DataType[]>
Equation_Computer::Equation_Computer::reorder_tensor_data( const DataType* orig_data, shared_ptr<vector<int>>  new_order_vec,
                                                           shared_ptr<vector<Index>> orig_index_blocks ) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "reorder_tensor_data" << endl;

  shared_ptr<vector<size_t>> rlen = get_sizes(orig_index_blocks);
  shared_ptr<vector<size_t>> new_order_st = make_shared<vector<size_t>>(new_order_vec->size());   
  size_t block_size = get_block_size(orig_index_blocks->begin(), orig_index_blocks->end());
  array<int,4> sort_options = {0,1,1,1};

  unique_ptr<DataType[]> reordered_data(new DataType[block_size]);

   for ( int ii = 0 ; ii != new_order_vec->size(); ii++) 
     new_order_st->at(ii) = new_order_vec->at(ii);

  Tensor_Sorter::Tensor_Sorter<double> TS ;

  size_t num_ids =  orig_index_blocks->size();
  if ( num_ids == 2) { 
    TS.sort_indices_2( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 3 ) {
    TS.sort_indices_3( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 4 ) {
    TS.sort_indices_4( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 5 ) {
    TS.sort_indices_5( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 6 ) {
    TS.sort_indices_6( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  } else if ( num_ids == 7 ) {
    TS.sort_indices_7( new_order_st, rlen, sort_options, orig_data, reordered_data.get() ) ;
  }

  return reordered_data;
}

#endif
