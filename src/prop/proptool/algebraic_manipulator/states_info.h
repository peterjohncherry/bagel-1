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
     std::shared_ptr<std::map< char ,int>> hole_range_map_;
     std::shared_ptr<std::map< char ,int>> elec_range_map_;

     long unsigned int elec_pnum_;
     long unsigned int hole_pnum_;

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
     std::shared_ptr<std::map< char, int>> hole_range_map() { return hole_range_map_; } 
     std::shared_ptr<std::map< char, int>> elec_range_map() { return elec_range_map_; } 

     long unsigned int elec_pnum() { return elec_pnum_ ; }  
     long unsigned int hole_pnum() { return hole_pnum_ ; }  
  
     // characeteristic numbers for determining  if < A | ..... | B > is zero without maps
     void set_elec_hole_pnums( std::shared_ptr<std::map<char, long unsigned int>> range_prime_map ) {

       hole_pnum_ = 1;
       for ( auto& elem : *hole_range_map_ ) 
         if ( elem.second != 0 )
           if ( elem.second < 5 ){ 
             hole_pnum_ *= (elem.second* range_prime_map->at(elem.first)); 
           } else { 
             hole_pnum_ *= pow(range_prime_map->at(elem.first), 5); 
           }

       elec_pnum_ = 1; 
       for ( auto& elem : *elec_range_map_ ) 
         if ( elem.second != 0 ) 
           if ( elem.second < 5 ){ 
             elec_pnum_ *= (elem.second* range_prime_map->at(elem.first)); 
           } else { 
             elec_pnum_ *= pow(range_prime_map->at(elem.first), 5); 
           }
     }

};

//Written strangely so can be compatible with states with multiple spin sectors.
//Would be better with base class and two derived classes.
template<class DataType>
class StatesInfo  { 
   
   public: 
     std::vector<int> target_state_nums_;
     std::map< std::string, std::shared_ptr<CIVecInfo<DataType>> > civec_info_map_;    
     std::map< int, std::shared_ptr<std::vector<std::string>> > state_civec_names_;

     bool multisector_; //TODO having civecs in multiple sectors should be dealt by templating class and having specialized function

     std::shared_ptr<std::map< char , long unsigned int >> range_prime_map_;

     StatesInfo(std::vector<int> Target_state_nums) : target_state_nums_(Target_state_nums), multisector_(false)  {}; // 
     ~StatesInfo(){}; 

     int spin_diff( std::string bra, std::string ket ) { return civec_info_map_.at(bra)->nalpha() - civec_info_map_.at(ket)->nalpha(); }         
  

     void add_state( const int nact, const int nele, const int state_number, std::shared_ptr<std::map<char, int>> elec_range_map,
                     std::shared_ptr<std::map<char, int>> hole_range_map ) {
        
       if ( state_civec_names_.find(state_number) == state_civec_names_.end() ){ 
         std::string civec_name = WickUtils::get_civec_name( state_number, elec_range_map->at('a') + hole_range_map->at('a'), elec_range_map->at('a'), elec_range_map->at('A'));
         civec_info_map_.emplace( civec_name, std::make_shared<CIVecInfo<DataType>>( state_number, civec_name, elec_range_map, hole_range_map ) ); 
         state_civec_names_.emplace(state_number, std::make_shared<std::vector<std::string>>(1, civec_name));
       } else { 
         std::cout << "state " << state_number << " already has it's civecs in the map !! " << std::endl; 
       }
     }
 
     //for spin-free only 
     std::shared_ptr<std::vector<std::string>> civec_names(int state_number) { return state_civec_names_.at(state_number); }

     int nalpha(int state_number) { return civec_info_map_.at(state_civec_names_.at(state_number)->front())->nalpha(); }
     int nbeta(int state_number) { return civec_info_map_.at(state_civec_names_.at(state_number)->front())->nbeta(); }
     int nact(int state_number) { return civec_info_map_.at(state_civec_names_.at(state_number)->front())->nact(); }
     int nele(int state_number) { return civec_info_map_.at(state_civec_names_.at(state_number)->front())->nele(); } 

     std::shared_ptr<std::vector< std::shared_ptr<CIVecInfo<DataType>>>> civec_info ( int state_number )  {
       return civec_info_map_.at(state_civec_names_.at(state_number));
     };
 
     std::shared_ptr<CIVecInfo<DataType>> civec_info ( std::string civec_name )  {
       return civec_info_map_.at(civec_name);
     };

     std::shared_ptr<std::map< char ,int>> hole_range_map( std::string civec_name )  {
       return civec_info_map_.at(civec_name)->hole_range_map();
     };
  
     std::shared_ptr<std::map< char ,int>> elec_range_map( std::string civec_name )  {
       return civec_info_map_.at(civec_name)->elec_range_map();
     };
  
};
#endif
