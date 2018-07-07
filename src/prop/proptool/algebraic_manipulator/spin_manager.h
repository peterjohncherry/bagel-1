#ifndef __SRC_PROP_PROPTOOL_Spin_Manager_H
#define __SRC_PROP_PROPTOOL_Spin_Manager_H

class Spin_Manager  {
  using pint_vec = std::vector<std::pair<int,int>>;
  using pstr_vec = std::vector<std::pair<std::string,std::string>>;
  
  public:

    int nact_el;
    int nact_orb;
    int nidxs;

    Spin_Manager(int nel, int norb, int num_idxs) :  nact_el(nel), nact_orb(norb),  nidxs(num_idxs) {}
    ~Spin_Manager(){};

    bool spinfree;
    int spin_max ;
    int spin_min ;

    std::shared_ptr<std::vector<std::pair<int,int>>> spin_sectors ;

    std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> >> spin_paths;
    
    std::shared_ptr<std::map< std::pair< int, std::pair<int,int> > , std::shared_ptr<std::vector<std::shared_ptr<std::vector<std::string>>>> >> spin_paths_str;

    void set_spin_info();
    
    std::shared_ptr<std::vector<std::pair<int,int>>> spin_sector_gen(int nact_el, int nact_orb);

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> all_spin_combinations(int length);

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> sector_constrained_spin_combs(int nalpha, int nbeta);

    bool sector_constrained_spin_combs(int nalpha, int nbeta , std::shared_ptr<std::vector<std::string>> gamma_comb);

    bool sector_constrained_spin_combs(int nalpha, int nbeta , std::shared_ptr<std::vector<bool>> gamma_comb);

    void bfunc(int xx, std::shared_ptr<std::vector<bool>> xvec, int nel, int length, std::shared_ptr<std::vector<std::shared_ptr<std::vector<bool>>>> all_sc);

    void get_string_spin_paths(); //should not be needed...

};
#endif
