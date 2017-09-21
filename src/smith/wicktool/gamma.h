 #include <src/smith/wicktool/WickUtils.h>
 #ifndef __SRC_SMITH_Gamma_H
 #define __SRC_SMITH_Gamma_H

// #include "WickUtils.h"
using namespace WickUtils;



class RDMderiv_new{ 
 
  public: 

  // variables
  bool spinfree;
  std::string name;

  std::shared_ptr<std::vector<bool>> full_aops ;
  std::shared_ptr<std::vector<std::string>> full_ids ;
  std::shared_ptr<std::vector<std::string>> full_id_ranges ;

  std::shared_ptr<std::vector<int>> ids_pos ;
  std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos; 

  std::shared_ptr<std::vector<int>> signs_all  ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> ids_pos_all ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<int,int>>>>> deltas_pos_all; 

  std::shared_ptr<std::map< char, int>>   op_order ;
  std::shared_ptr<std::map< std::string, int>>  idx_order ;
  
  int  my_sign ;

  RDMderiv_new(){};
  RDMderiv_new( std::shared_ptr<std::vector<bool>> full_aops_in , std::shared_ptr<std::vector<std::string>> full_ids_in ,  std::shared_ptr<std::vector<std::string>> full_id_ranges_in ,
                std::shared_ptr<std::vector<int>> ids_pos_in ,  std::shared_ptr<std::vector<std::pair<int,int>>> deltas_pos_in, int sign_in ) :
                full_aops(full_aops_in), full_ids(full_ids_in), full_id_ranges(full_id_ranges_in), ids_pos(ids_pos_in), deltas_pos(deltas_pos_in), my_sign(sign_in) {};

  ~RDMderiv_new(){};

 // functions
  void initialize(std::shared_ptr<std::vector<bool>> ac_init,std::shared_ptr<std::vector<std::string>> ids_init,
                  std::shared_ptr<std::vector<std::string>> id_ranges);

  void initialize( std::shared_ptr<std::vector<bool>> ac_init,  std::shared_ptr<std::vector<std::string>> id_ranges, std::shared_ptr<std::vector<std::string>> ids_init,
                   std::shared_ptr<std::vector<std::pair<int,int>>> deltas_init, int sign = 1);


 

  void swap(std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<pint_vec> deltas_pos, int ii, int jj, int kk );

  void generic_reordering(std::shared_ptr<std::vector<int>> new_order );

  std::shared_ptr<std::vector<std::pair<int,int>>> Standardize_delta_ordering(std::shared_ptr<std::vector<std::pair<int,int>>> delta_ids ) ;

  void norm_order();
  
  void alt_order();

  void norm_order_recursive(std::shared_ptr<std::vector<std::shared_ptr<RDMderiv_new>>> rdm_vec );
 
  void swap_recursive(std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<pint_vec> deltas_pos, int ii, int jj, int kk,
                      std::shared_ptr<std::vector<std::shared_ptr<RDMderiv_new>>> rdm_vec  );

  bool ordering(std::pair<std::string,std::string> a, std::pair<std::string,std::string> b) {
       return( (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false ); };

    
  bool gamma_survives(std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<std::vector<std::string>> id_ranges);
  
  bool gamma_constraints(std::shared_ptr<std::vector<int>> ids_pos, std::shared_ptr<std::vector<std::string>> id_ranges, std::shared_ptr<std::vector<bool>> full_aops);
  

  std::string get_Aname(std::shared_ptr<std::vector<std::string>> full_idxs, std::shared_ptr<std::vector<std::string>> full_idx_ranges,
                        std::shared_ptr<std::vector<std::pair<int,int>>> all_ctrs_pos );

  std::string get_gamma_name(std::shared_ptr<std::vector<std::string>> full_idx_ranges,  std::shared_ptr<std::vector<bool>> orig_aops,
                             std::shared_ptr<std::vector<int>> all_idxs_pos );


};




class RDMderiv{ 
 
  public: 
    RDMderiv(){};
    ~RDMderiv(){};

 // functions
  void initialize( std::shared_ptr<std::vector<bool>> ac_init, std::shared_ptr<std::vector<std::string>> ids_init,
                   std::shared_ptr<std::vector<std::pair<std::string,std::string>>> deltas_init, int sign = 1);
 
  void initialize(std::shared_ptr<std::vector<bool>> ac_init, std::shared_ptr<std::vector<std::string>> ids_init, std::shared_ptr<std::vector<std::string>> ids_spins);

  void swap(std::shared_ptr<std::vector<bool>> a_vector, std::shared_ptr<std::vector<std::string>> ids,
            std::shared_ptr<std::vector<std::pair<std::string, std::string>>> dlist, int ii, int jj, int kk  );
 
  void norm_order();
 
  //this routine can be used for arbitrary ordering; for testing only
  void generic_reordering(std::shared_ptr<std::vector<bool>> new_order );

  std::shared_ptr<std::vector<std::pair<std::string,std::string>>> Standardize_delta_ordering(std::shared_ptr<std::vector<std::pair<std::string,std::string>>> delta_ids ) ;
 

  // variables
  std::shared_ptr<std::vector<std::string>> orig_ids ;
  std::shared_ptr<std::vector<bool>> aops ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>>   allops ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> allids ;
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string,std::string>>>>> alldeltas; 
  std::shared_ptr<std::vector<int>> allsigns  ;
  bool spinfree;
  std::shared_ptr<std::map< char, int>>   op_order ;
  std::shared_ptr<std::map< std::string, int>>  idx_order ;
  std::shared_ptr<std::map< std::string, bool>> aop_map   ;
  
  bool ordering(std::pair<std::string,std::string> a, std::pair<std::string,std::string> b) { bool tim = (idx_order->at(a.first) < idx_order->at(b.first) ) ? true : false; return tim; };
};

//would be called gamma, but apparently that's part of stl
class alt_RDMderiv : public RDMderiv { 
 
  public: 
    alt_RDMderiv(){};
    ~alt_RDMderiv(){};

 // functions
  void alt_order();
  
  std::shared_ptr<pint_vec> get_spinsector_path(std::shared_ptr<pstr_vec> gamma_spin_ids , std::pair<int,int> spinsector ) ;
  std::shared_ptr<std::vector<std::string>> get_ranges(std::shared_ptr<std::vector<std::string>> spin_ids ) ;


  //note that this will remove the spins from the end of the delta_spin_ids, whereas get_spin (defined above) just returns the spins without editing spin_ids
  std::shared_ptr<std::vector<std::pair<std::string,std::string>>> strip_delta_spins(std::shared_ptr<std::vector<std::pair<std::string,std::string>>> deltas_spin_ids ) ;


  // variables
  std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::pair<std::string,std::string>>>>> alldeltas_spins; 


};

class Single_Gamma {

   public:
     Single_Gamma(){};
     ~Single_Gamma(){};

   //variables
   std::shared_ptr<pint_vec> spinsector_path;
   std::shared_ptr<std::vector<std::string>> gamma_spin_ids;
   bool calculated;  

   //functions
   std::shared_ptr<double> compute(); 
   std::shared_ptr<double> data;
}; 
#endif
