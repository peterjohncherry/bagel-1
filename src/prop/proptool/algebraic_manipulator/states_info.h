#ifndef __SRC_PROP_PROPTOOL_Target_Info_H
#define __SRC_PROP_PROPTOOL_Target_Info_H

#include <src/prop/proptool/proputils.h>
template<class DataType>
class CIVecInfo  {
   public:
     CIVecInfo(){}
     ~CIVecInfo(){};
};


template<> class CIVecInfo<double> {
 
   protected: 
     const int nalpha_;
     const int nbeta_;
     const int nact_;
     const int state_number_;
     const std::string name_;

   public:

     bool sparse_;

     CIVecInfo(int nalpha, int nbeta , int nact, int state_number, std::string name  ) :
     nalpha_(nalpha), nbeta_(nbeta), nact_(nact), state_number_(state_number), name_(name), sparse_(false) {};
  
     ~CIVecInfo(){};
      
     int nele()   { return nalpha_+nbeta_; };
     int nalpha() { return nalpha_; };
     int nbeta()  { return nbeta_;  };
     int nact()   { return nact_;   };
     std::string name() { return name_;};
     int state_num() { return state_number_;};

};

template<> class CIVecInfo<std::complex<double>>  {
 
   private : 
     const int nact_;
     const int nele_;
     const int state_number_;
     std::map< std::string, std::shared_ptr<CIVecInfo<double>> > civec_info_map;    

   public :
     CIVecInfo(int nact, int nele, int state_number) :
     nact_(nact), nele_(nele), state_number_(state_number)  {

       int min_alpha = ( nele_ <= nact_ ) ? 0 : nele_ - nact_ ;
       int max_alpha = ( nact_ <= nele_ ) ? nact_ : nele_ ;
       for ( int nalpha = min_alpha ;  nalpha != max_alpha; nalpha++ ) {

         std::string civec_name = WickUtils::get_civec_name( state_number, nact_, nalpha, nele_ - nalpha);
         civec_info_map.emplace( civec_name, std::make_shared<CIVecInfo<double>>(nalpha, nele_ - nalpha , nact_, state_number_, civec_name )); 

       }
     };
  
     ~CIVecInfo(){};

     int nele() { return nele_; }
     int nact() { return nact_;   }
     int state_number() { return nact_; }

};
//Written strangely so can be compatible with states with multiple spin sectors.
//Would be better with base class and two derived classes.
template<class DataType>
class StatesInfo  { 
   
   public: 
     std::vector<int> target_state_nums_;
     std::map< std::string, std::shared_ptr<CIVecInfo<double>> > civec_info_map;    

     //for spin-free only
     std::map< int, std::string > state_name;

     StatesInfo(std::vector<int> Target_state_nums) : target_state_nums_(Target_state_nums) {};
     ~StatesInfo(){}; 

     //should go to alternate function in relativistic case
     void add_state( const int nalpha, const int nbeta , const int nact, const int state_number )  {
         
        std::string civec_name = WickUtils::get_civec_name( state_number, nact, nalpha, nbeta);
        civec_info_map.emplace( civec_name, std::make_shared<CIVecInfo<double>>(nalpha, nbeta , nact, state_number, civec_name )); 
        state_name.emplace( state_number, civec_name ); 

     }        

     void add_state( const int nact, const int nele, const int state_number ) {

       int min_alpha = ( nele <= nact ) ? 0 : nele - nact ;
       int max_alpha = ( nact <= nele ) ? nact : nele ;

       for ( int nalpha = min_alpha ;  nalpha != max_alpha; nalpha++ ) {
         std::string civec_name = WickUtils::get_civec_name( state_number, nact, nalpha, nele - nalpha);
         civec_info_map.emplace( civec_name, std::make_shared<CIVecInfo<double>>(nalpha, nele - nalpha , nact, state_number, civec_name )); 
       }
     }
     
     int spin_diff( std::string bra, std::string ket ) { return civec_info_map.at(bra)->nalpha() - civec_info_map.at(ket)->nalpha(); }         
 
     //for spin-free only 
     std::string name(int state_number)   { return state_name.at(state_number); }

     int nalpha(int state_number) { return civec_info_map.at(state_name.at(state_number))->nalpha(); }
     int nbeta(int state_number) { return civec_info_map.at(state_name.at(state_number))->nbeta(); }
     int nact(int state_number) { return civec_info_map.at(state_name.at(state_number))->nact(); }
     int nele(int state_number) { return civec_info_map.at(state_name.at(state_number))->nele(); } 

     std::shared_ptr<CIVecInfo<double>> civec_info ( int state_number )  { return civec_info_map.at(state_name.at(state_number)); };
  
};
#endif
