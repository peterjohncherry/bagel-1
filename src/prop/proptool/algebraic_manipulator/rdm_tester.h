#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

class RDM_Tester {

  public:  
    
    std::shared_ptr<StatesInfo_Base> states_info_;
    std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo_Base>>>& gamma_map_;

    RDM_Tester( std::shared_ptr<StatesInfo_Base> states_info,  std::shared_ptr<std::map<std::string, std::shared_ptr<GammaInfo_Base>>>& gamma_map ) :
                states_info_(states_info), gamma_map_(gamma_map) {}  
   ~RDM_Tester(){}
 
    void make_gamma_info( int order, int bra_name, int ket_name  );
 
    void make_gamma_info( int order, std::string bra_name, std::string ket_name  );
};
