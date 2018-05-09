#include <bagel_config.h>
#include <src/prop/proptool/algebraic_manipulator/op_info.h>

using namespace std;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MultiOp_Info::MultiOp_Info( vector<string>& op_list, vector<char>& op_trans_list, shared_ptr<vector<vector<int>>> op_state_ids ) { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "MultiOp_Info::MultiOp_Info" << endl;

  op_order_ = make_shared<vector<int>>(op_list.size());
  canonical_order_ = true;
 
  if ( op_order_->size() > 1 ) {
    iota( op_order_->begin(), op_order_->end(), 0);
    sort( op_order_->begin(), op_order_->end(), [ &op_list ] ( int i1, int i2) { return (bool)( op_list[i1] < op_list[i2]); });  
  
    for ( vector<int>::iterator oo_it = op_order_->begin()+1; oo_it !=op_order_->end(); oo_it++ ) {
      if ( *oo_it != *(oo_it-1)+1 ){
        canonical_order_ = false;
        break;
      }
    }
  }

  op_name_       = "";
  op_state_name_ = "";
  op_full_name_  = "";

  shared_ptr<vector<shared_ptr<vector<int>>>> state_id_list =  make_shared<vector<shared_ptr<vector<int>>>>();

  op_info_vec_ = make_shared<vector<shared_ptr<Op_Info>>>( op_list.size() );
  vector<shared_ptr<Op_Info>>::iterator oiv_it = op_info_vec_->begin();

  for ( int ii = 0 ; ii != op_list.size(); ii++, oiv_it++ ) {

    string single_op_full_name = op_list[ii] ; // Because cannot be set equal to char on initialization...
    string single_op_name = single_op_full_name;

    if (op_state_ids->at(ii).size() > 0 ) {
      single_op_full_name +=  "_{"; 
      for( int jj = 0; jj != (*op_state_ids)[ii].size(); jj++ ) {
        single_op_full_name += to_string((*op_state_ids)[ii][jj]); 
      }
      single_op_full_name += "}"; 
    }

    string single_op_state_name = single_op_full_name;

   if (op_trans_list.size() > 0 ) {
      single_op_full_name +=  "^{"; 
      single_op_full_name += op_trans_list[ii]; 
      single_op_full_name += "}"; 
    }
    *oiv_it = make_shared<Op_Info>( single_op_name, single_op_state_name, single_op_full_name, make_shared<vector<int>>( (*op_state_ids)[ii] ) , op_trans_list[ii] ); 

    op_name_       += single_op_name;
    op_state_name_ += single_op_state_name;
    op_full_name_  += single_op_full_name;

  }

  num_ops_ = op_info_vec_->size();  
  state_ids_list_ = make_shared<vector<shared_ptr<vector<int>>>>(num_ops_);    

  transformations_ = make_shared<vector<char>>(num_ops_);    
  vector<char>::iterator t_it = transformations_->begin(); 
  vector<shared_ptr<vector<int>>>::iterator sil_it = state_ids_list_->begin(); 
  for ( vector<shared_ptr<Op_Info>>::iterator oi_it = op_info_vec_->begin(); oi_it != op_info_vec_->end(); oi_it++ , t_it++, sil_it++ ){
    *t_it = (*oi_it)->transformation_; 
    *sil_it = (*oi_it)->state_ids_; 
  }

  if ( !canonical_order_ ) { 
    vector<string> op_list_canonical(num_ops_);
    vector<char> op_trans_list_canonical(num_ops_);
    shared_ptr<vector<vector<int>>> state_ids_canonical = make_shared<vector<vector<int>>>(num_ops_);

    vector<string>::iterator      olc_it  = op_list_canonical.begin();
    vector<char>::iterator        otlc_it = op_trans_list_canonical.begin();
    vector<vector<int>>::iterator sic_it  = state_ids_canonical->begin();
    for ( int pos : *op_order_ ) {
      *olc_it++ = op_list[pos];
      *otlc_it++ = op_trans_list[pos];
      *sic_it++ = (*op_state_ids)[pos];
    } 
    op_info_canonical_ = make_shared<MultiOp_Info>( op_list_canonical, op_trans_list_canonical, state_ids_canonical ); 
  }   
}
