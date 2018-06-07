#include <bagel_config.h>
#ifdef COMPILE_SMITH
#include <src/prop/proptool/proputils.h>
 // #include "wickutils.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Because always want to keep real and imaginary factors seperate
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void WickUtils::pair_fac_mult( const std::pair<double,double>& factor_fixed , std::pair<double,double>& factor_changing  ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::pair_fac_mult" << endl;

  double re_buff = factor_changing.first;
  double im_buff = factor_changing.second; // don't need two buffers, but for clarities sake...

  factor_changing.first = re_buff*factor_fixed.first - im_buff*factor_fixed.second;
  factor_changing.second = im_buff*factor_fixed.first + re_buff*factor_fixed.second;
 
  return; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//same as  fvec_cycle, but allows skipping. Should be included everywhere to guard against max==min problem};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper(shared_ptr<vector<int>> forvec, shared_ptr<vector<int>> max , shared_ptr<vector<int>> min ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec->size()-1; ii!=-1 ; ii--) {
    if ( max->at(ii) == min->at(ii) ) {
      if ( ii == 0 )
        return false;
    } else if (forvec->at(ii) == max->at(ii)) {
      if (ii == 0) 
        return false;    
      forvec->at(ii) = min->at(ii);
    } else {
      forvec->at(ii) = forvec->at(ii)+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//NO SHARED PTR VERSION, SHOULD BE BETTER THAN SHARED POINTERS
//same as fvec_cycle, but allows skipping. Should be included everywhere to guard against max==min problem};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper( vector<int>& forvec, vector<int>& max, vector<int>& min ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = forvec.size()-1; ii!=-1 ; ii--) {
    if ( max[ii] == min[ii] ) {
      if ( ii == 0 )
        return false;
    } else if (forvec[ii] == max[ii]) {
      if (ii == 0) 
        return false;    
      forvec[ii] = min[ii];
    } else {
      forvec[ii] = forvec[ii]+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
//REVERSE ITERATOR VERSION, SHOULD BE FASTEST, BUT CHECK 
//TODO TEMPLATE (AND MAYBE MODIFY) SO IT WORKS FOR DIFFERENT TYPES OF CONTAINER
//same as fvec_cycle, but allows skipping. Should be included everywhere to guard against max==min problem};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper( vector<int>& forvec, vector<int>::reverse_iterator max_it, vector<int>::reverse_iterator min_it ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(vector<int>::reverse_iterator f_it = forvec.rbegin(); f_it !=forvec.rend(); f_it++,  max_it++, min_it++ ) {
    if ( *max_it == *min_it ) {
      if ( f_it == (forvec.rend()-1) )
        return false;
    } else if (*f_it == *max_it ) {
      if ( f_it == (forvec.rend()-1) )
        return false;    
      *f_it = *min_it;
    } else {
      (*f_it) += 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::fvec_cycle_skipper_f2b( vector<int>& forvec, vector<int>& max , vector<int>& min ) {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int ii = 0; ii!=forvec.size(); ii++) {
    if ( max[ii] == min[ii] ) {
      if ( ii == forvec.size()-1 )
        return false;
    } else if (forvec[ii] == max[ii]) {
      if (ii == forvec.size()-1) 
        return false;    
      forvec[ii] = min[ii];
    } else {
      forvec[ii] = forvec[ii]+ 1;
      break;
    }
  }
  return true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> WickUtils::reorder_vector(vector<int>& neworder , const vector<int>& origvec ) {
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  vector<int> newvec(origvec.size());
  vector<int>::iterator newvec_it = newvec.begin();

  for( int pos : neworder )
     *newvec_it++ = origvec[pos];

  return make_shared<vector<int>>(newvec);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_civec_name(const int state_num, const int norb, const int nalpha, const int nbeta)  { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = to_string(state_num) + "_["+ to_string(norb)+"o{" + to_string(nalpha) + "a," + to_string(nbeta) + "b}]" ;
  return name ;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_gamma_name( const vector<string>& full_idx_ranges,  const vector<bool>& aops_vec,
                                  const vector<int>& idxs_pos, string Bra_name, string Ket_name ){
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "WickUtils::get_gamma_name" << endl; 
  string  name;
 
  if (idxs_pos.size() == 0 ) {
     name = "ID" ;
  } else {  
    name = "<" + Bra_name + "|_(";
    for (int pos : idxs_pos ) 
      name+=full_idx_ranges[pos][0];
    
    name+='_';
    for (int pos : idxs_pos ) {
      if(aops_vec[pos]){ 
        name += '1';
      } else {
        name += '0';
      }
    } 
    name += ")_|" + Ket_name + ">";
  }
  
  return name;
};
////////////////////////////////////////////////////////////////////////////////////////////////////
// Check if the range block will survive without contraction with another
// range block (i.e. see if as many particles are being created in a given
// range as are destroyed) 
////////////////////////////////////////////////////////////////////////////////////////////////////
bool WickUtils::RangeCheck(const vector<string>& id_ranges, const vector<bool>& aops ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::RangeCheck" << endl;

  vector<string> diff_rngs(1, id_ranges[0] );
  vector<int> updown(1, ( aops[0] ? 1 : -1)  );

  for ( int jj = 1;  jj !=id_ranges.size() ; jj++){
    int ii = 0;

    string rng = id_ranges[jj];

    do {
      if(rng == diff_rngs[ii]){
        if ( aops[jj]){
          updown[ii]+=1;
        } else {         
          updown[ii]-=1;
        }
        break;
      }
      if ( ii == diff_rngs.size()-1){
        diff_rngs.push_back(rng);
        if ( aops[jj]){
          updown.push_back(1);  
        } else {
          updown.push_back(-1); 

        }
        break;
      }

      ii++;
    } while (true);
  } 

  for (int ac : updown ) 
    if (ac != 0 )
      return false;
 
  return true;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
// get range for ctrtensorpart  
////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_ctp_name( const string op_state_name, const vector<string>& idxs,
                                const vector<string>& id_ranges, const vector<pair<int,int>>& ctrs_pos ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//  cout << "WickUtils::get_ctp_name" << endl; 

  string ctp_name = op_state_name+"_";

  for (string id : id_ranges ) 
    ctp_name += id[0];

  if (ctrs_pos.size() != 0 ){ 
    ctp_name += "_";
    vector<pair<int,int>> ctrs_buff =  ctrs_pos;
    shared_ptr<vector<pair<int,int>>> ctrs_standard = standardize_delta_ordering_generic(ctrs_buff, idxs); 
    for ( pair<int,int>& ctr : *ctrs_standard )  
      ctp_name += to_string(ctr.first)+to_string(ctr.second);
  }
 
  return ctp_name;

}
////////////////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<pair<int,int>>>
WickUtils::standardize_delta_ordering_generic( const vector<pair<int,int>>& deltas_pos, const vector<string>& idxs  ) {
////////////////////////////////////////////////////////////////////////////////////////////////////
//cout << "GammaGenerator::Standardize_delta_ordering_generic" << endl;
//TODO must order by indexes, not just by initial position
//     If one of the indexes is X, cannot "just" not contract
//     must also account for reordering ; T_{ijkl} = ... + <I| ijklmnop | J> A_{mnop} delta_{lm}
//     T_{ijkl} = .... + < I | ijklmnop | J> A_{lmop} 
//     so must flip order of A based on T (X) position

 shared_ptr<vector<pair<int,int>>> new_deltas_pos;

 if (deltas_pos.size() > 1 ) {
   new_deltas_pos =  make_shared <vector<pair<int,int>>>( deltas_pos.size());
   vector<int> posvec(deltas_pos.size(),0);
   for (int ii = 0 ; ii != deltas_pos.size() ; ii++)
     for (int jj = 0 ; jj != deltas_pos.size() ; jj++)
       if (idxs[deltas_pos[ii].first] > idxs[deltas_pos[jj].first] )
         posvec[ii]++;

   for (int ii = 0 ; ii != deltas_pos.size() ; ii++)
     new_deltas_pos->at(posvec[ii]) = deltas_pos[ii];

 } else {
   vector<pair<int,int>> new_deltas_pos_raw(deltas_pos);
   new_deltas_pos = make_shared<vector<pair<int,int>>>( new_deltas_pos_raw ) ;

 }
  return new_deltas_pos;
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::range_to_prime(char range ) {
//////////////////////////////////////////////////////////////////////////////////
  switch (range){
    case 'v' : return 2;
    case 'V' : return 3;
    case 'c' : return 5;
    case 'C' : return 7;
    case 'a' : return 11;
    case 'A' : return 13;
    default : 
      throw std::logic_error( " unknown range " ); return 9999999;
  }
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::range_to_prime_spinfree(char range ) {
//////////////////////////////////////////////////////////////////////////////////
  switch (range){
    case 'v' : return 2;
    case 'c' : return 3;
    case 'a' : return 5;
    default : 
      throw std::logic_error( " unknown range " ); return 9999999;
  }
}
//////////////////////////////////////////////////////////////////////////////////
unsigned int WickUtils::get_block_hash( const std::vector<std::string>& block  ) {
//////////////////////////////////////////////////////////////////////////////////
 
  int max_range_num = 13; // should set this in main proptool
  int range_id = 0;
  int ii = 0;
  for ( vector<string>::const_reverse_iterator b_it = block.crbegin() ; b_it != block.crend() ;  b_it++, ii++ ) 
     range_id += range_to_prime( (*b_it)[0] ) * pow( max_range_num, ii ); 
 
  assert( range_id >  -1 );

  return range_id; 
}
//////////////////////////////////////////////////////////////////////////////////////
vector<char> WickUtils::strvec_to_chrvec( const vector<string>& strvec ) {
//////////////////////////////////////////////////////////////////////////////////////
 
   vector<char> chrvec(strvec.size());
   vector<char>::iterator cv_it = chrvec.begin();
   for( vector<string>::const_iterator sv_it = strvec.begin(); sv_it != strvec.end(); sv_it++, cv_it++ )
     *cv_it = (*sv_it)[0];
  
   return chrvec  ; 
} 
//////////////////////////////////////////////////////////////////////////////////////
vector<string> WickUtils::chrvec_to_strvec( const vector<char>& chrvec ) {
//////////////////////////////////////////////////////////////////////////////////////
 
  vector<string> strvec(chrvec.size());
  vector<char>::const_iterator cv_it = chrvec.begin();
  for ( vector<string>::iterator sv_it = strvec.begin(); sv_it != strvec.end(); sv_it++, cv_it++ ) {
    string s = "";
    s += *cv_it;
    *sv_it = s;
  }

  return strvec;
} 
//////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<int>> WickUtils::get_ascending_order( const vector<int>& scrambled_vec ) {
//////////////////////////////////////////////////////////////////////////////////////
 
  shared_ptr<vector<int>> new_order = make_shared<vector<int>>(scrambled_vec.size());
  iota(new_order->begin(), new_order->end() , 0 );
  sort ( new_order->begin(), new_order->end(), [&scrambled_vec]( int& i1, int& i2 ){ return (bool)( scrambled_vec[i1] < scrambled_vec[i2] );}); 

  return new_order;
} 
#endif
