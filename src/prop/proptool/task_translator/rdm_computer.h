#ifndef __SRC_PROP_PROPTOOL_RDM_COMPUTER_H
#define __SRC_PROP_PROPTOOL_RDM_COMPUTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_sorter.h>
#include <src/prop/proptool/tensor_and_ci_lib/tensor_arithmetic.h>
#include <src/prop/proptool/algebraic_manipulator/gamma_info.h>

namespace bagel {

namespace RDM_Computer { 
class RDM_Computer { 
 
    public: 

      std::unique_ptr<RDM_Tester> rdm_tester_;
      std::shared_ptr<std::map< std::string , std::shared_ptr<SMITH::Tensor_<double>> >> civec_data_map_;
      std::shared_ptr<std::map< std::string , std::shared_ptr<SMITH::Tensor_<double>> >> gamma_data_map_;

      RDM_Computer( std::shared_ptr<std::map< std::string , std::shared_ptr<SMITH::Tensor_<double>> >> civec_data_map,
                    std::shared_ptr<std::map< std::string , std::shared_ptr<SMITH::Tensor_<double>> >> gamma_data_map ):
                    civec_data_map_(civec_data_map), gamma_data_map_(gamma_data_map) {}
      ~RDM_Computer(){};
      
};

}
}
#endif
