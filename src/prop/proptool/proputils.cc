#include <bagel_config.h>
#include <src/prop/proptool/proputils.h>
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string WickUtils::get_civec_name(const int state_num, const int norb, const int nalpha, const int nbeta)  { 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  string name = to_string(state_num) + "_["+ to_string(norb)+"o{" + to_string(nalpha) + "a," + to_string(nbeta) + "b}]" ;
  return name ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string WickUtils::get_det_name( char name_range1, int num_range1, char name_range2, int num_range2, int norb ){
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   char num_range1_char = '0'+ num_range1;
   char num_range2_char = '0'+ num_range2;
   char norb_char = '0'+ norb;

   string det_name = "(";
   det_name += num_range1_char;
   det_name += name_range1;
   det_name += ", ";
   det_name += num_range2_char;
   det_name += name_range2;
   det_name += " | ";
   det_name +=  norb_char;
   det_name += "o ) ";  

   return det_name;
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
std::vector<int> WickUtils::get_ascending_order( const std::vector<int>& scrambled_vec ) {
//////////////////////////////////////////////////////////////////////////////////////
    std::vector<int> new_order(scrambled_vec.size());
    std::iota( new_order.begin(), new_order.end() , 0 );
    std::sort( new_order.begin(), new_order.end(), [&scrambled_vec]( int& i1, int& i2 ){ return (bool)( scrambled_vec[i1] < scrambled_vec[i2] );}); 
    return new_order;
}; 

