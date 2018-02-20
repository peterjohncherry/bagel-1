#ifndef __SRC_PROP_PROPTOOL_Target_Info_H
#define __SRC_PROP_PROPTOOL_Target_Info_H

#include <src/prop/proptool/proputils.h>

// for ranges, lower case indicates an alpha range, and upper case indicates a beta range
// the letter specifies the spatial range.
template< typename DataType >
class CIVecInfo {
 
   protected: 
     const int nalpha_;
     const int nbeta_;
     const int nact_;
     const int nele_;
     const int state_number_;
     const std::string name_;
     std::shared_ptr<std::map< char ,int>>  hole_range_map_;
     std::shared_ptr<std::map< char ,int>>  elec_range_map_;

   public:

     bool sparse_;

     CIVecInfo( int state_number, std::string name , std::shared_ptr<std::map< char , int >> elec_range_map,
                std::shared_ptr<std::map< char , int >> hole_range_map ) :
     nalpha_( elec_range_map->at('a') ), nbeta_( elec_range_map->at('A') ),
     nact_(nalpha_ + hole_range_map->at('a')), nele_(nalpha_ + elec_range_map->at('A')), 
     state_number_(state_number), name_(name), elec_range_map_(elec_range_map),
     hole_range_map_(hole_range_map), sparse_(false) {};

     ~CIVecInfo(){};
      
     int nele()   { return nalpha_+nbeta_; };
     int nalpha() { return nalpha_; };
     int nbeta()  { return nbeta_;  };
     int nact()   { return nact_;   };
     std::string name() { return name_;};
     int state_num() { return state_number_;};
     bool sparse() { return sparse_;}
     std::shared_ptr<std::map< char ,int>> hole_range_map() { return hole_range_map_; } 
     std::shared_ptr<std::map< char ,int>> elec_range_map() { return elec_range_map_; } 

};

//Written strangely so can be compatible with states with multiple spin sectors.
//Would be better with base class and two derived classes.
template<class DataType>
class StatesInfo  { 
   
   public: 
     std::vector<int> target_state_nums_;
     std::map< std::string, std::shared_ptr<CIVecInfo<DataType>> > civec_info_map;    

     //for spin-free only
     std::map< int, std::string > state_name_;

     StatesInfo(std::vector<int> Target_state_nums) : target_state_nums_(Target_state_nums) {};
     ~StatesInfo(){}; 

     int spin_diff( std::string bra, std::string ket ) { return civec_info_map.at(bra)->nalpha() - civec_info_map.at(ket)->nalpha(); }         
 
     void add_state( const int nact, const int nele, const int state_number, std::shared_ptr<std::map<char, int>> elec_range_map,
                     std::shared_ptr<std::map<char, int>> hole_range_map ) {
       std::string civec_name = WickUtils::get_civec_name( state_number, elec_range_map->at('a') + hole_range_map->at('a'), elec_range_map->at('a'), elec_range_map->at('A'));
       state_name_.emplace(state_number, civec_name);
       civec_info_map.emplace( civec_name, std::make_shared<CIVecInfo<DataType>>( state_number, civec_name, elec_range_map, hole_range_map ) ); 
     }
 
     //for spin-free only 
     std::string name(int state_number)   { return state_name_.at(state_number); }

     int nalpha(int state_number) { return civec_info_map.at(state_name_.at(state_number))->nalpha(); }
     int nbeta(int state_number) { return civec_info_map.at(state_name_.at(state_number))->nbeta(); }
     int nact(int state_number) { return civec_info_map.at(state_name_.at(state_number))->nact(); }
     int nele(int state_number) { return civec_info_map.at(state_name_.at(state_number))->nele(); } 

     std::shared_ptr<CIVecInfo<DataType>> civec_info ( int state_number )  { return civec_info_map.at(state_name_.at(state_number)); };
  
};
#endif
