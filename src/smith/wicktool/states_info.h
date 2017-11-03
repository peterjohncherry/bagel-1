#ifndef __SRC_SMITH_Target_Info_H
#define __SRC_SMITH_Target_Info_H

template<class DataType>
class StatesInfo  { 
  
   public: 
     std::map<int,int> num_alpha_map;    
     std::map<int,int> num_beta_map;    
     std::map<int,int> num_orb_map;    
     std::vector<int> target_states; 
     
     StatesInfo (std::shared_ptr<std::vector<int>> target_states_in ) : target_states(*target_states_in) {};

     int nalpha(int state_num) { return num_alpha_map.at(state_num); }
     int nbeta(int state_num)  { return num_beta_map.at(state_num); }
     int norb(int state_num)   { return num_orb_map.at(state_num); }

     //should always be the same...
     int nele(int state_num)   { return ( num_alpha_map.at(state_num) + num_beta_map.at(state_num) ); }
      

};
#endif
