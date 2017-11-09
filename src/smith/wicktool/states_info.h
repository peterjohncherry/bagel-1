#ifndef __SRC_SMITH_Target_Info_H
#define __SRC_SMITH_Target_Info_H

#include <src/smith/wicktool/WickUtils.h>

template<class DataType>
class CIVecInfo  { 
  
   public: 
     const int nalpha_;    
     const int nbeta_;
     const int nact_;    
     const int state_number_;    
     const std::string name_; 

     CIVecInfo(int nalpha, int nbeta , int nact, int state_number, std::string name  ) :
     nalpha_(nalpha), nbeta_(nbeta), nact_(nact), state_number_(state_number), name_(name) {}; 
  
     ~CIVecInfo(){};

     int nele()   { return nalpha_+nbeta_; };
     int nalpha() { return nalpha_; };
     int nbeta()  { return nbeta_;  };
     int nact()   { return nact_;   };
     std::string name() { return name_;};
 
     //spin-free only
     int state_num() { return state_number_;};

};

//Written strangely so can be compatible with states with multiple spin sectors.
//Would be better with base class and two derived classes.
template<class DataType>
class StatesInfo  { 
  
   public: 
     std::vector<int> target_state_nums_;
     std::map< std::string, std::shared_ptr<CIVecInfo<DataType>> > civec_info_map;    

     //for spin-free only
     std::map< int, std::string > nr_state_name;

     StatesInfo(std::vector<int> Target_state_nums) : target_state_nums_(Target_state_nums) {};
     ~StatesInfo(){}; 

     //should go to alternate function in relativistic case
     void add_state( const int nalpha, const int nbeta , const int nact, const int state_number )  {
          
        std::string civec_name = WickUtils::get_civec_name( state_number, nact, nalpha, nbeta);

        nr_state_name.emplace( state_number, civec_name ); 
        civec_info_map.emplace( civec_name, std::make_shared<CIVecInfo<DataType>>(nalpha, nbeta , nact, state_number, civec_name )); 

     };

     //for spin-free only 
     std::string name(int state_number)   { return nr_state_name.at(state_number); }

     int nalpha(int state_number) { return civec_info_map.at(nr_state_name.at(state_number))->nalpha(); }
     int nbeta(int state_number)  { return civec_info_map.at(nr_state_name.at(state_number))->nbeta(); }
     int nact(int state_number)   { return civec_info_map.at(nr_state_name.at(state_number))->nact(); }
     int nele(int state_number)   { return civec_info_map.at(nr_state_name.at(state_number))->nele(); } 

     std::shared_ptr<CIVecInfo<DataType>> civec_info ( int state_number )  { return civec_info_map.at(nr_state_name.at(state_number)); };
  
};
#endif
