#ifndef __SRC_PROP_PROPTOOL_TENS_AND_CILIB_CI_TYPE_CONVERTER_H
#define __SRC_PROP_PROPTOOL_TENS_AND_CILIB_CI_TYPE_CONVERTER_H

#include <src/smith/tensor.h>
#include <src/smith/indexrange.h>
#include <src/ci/fci/civec.h>
#include <src/util/f77.h>

namespace bagel {

namespace CI_Type_Converter { 
template<typename DataType>
class CI_Type_Converter { 
  public: 

    // TODO For different conversions modify these using statements (why not typedef?)
    using Index = SMITH::Index;
    using IndexRange = SMITH::IndexRange;
    using Tensor = typename std::conditional<std::is_same<DataType,double>::value,
                                             SMITH::Tensor_<double>,SMITH::Tensor_<std::complex<double>>>::type;

    std::shared_ptr<std::map< std::string, std::shared_ptr<Tensor>>> civec_data_map_;
    std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map_; 
    std::shared_ptr<std::map< std::string, std::shared_ptr<const Determinants>>> bagel_determinant_map_;


    CI_Type_Converter(){};
    CI_Type_Converter( std::shared_ptr<std::map< std::string, std::shared_ptr<IndexRange>>> range_conversion_map ): 
                       range_conversion_map_(range_conversion_map),
                       civec_data_map_( std::make_shared<std::map< std::string, std::shared_ptr<Tensor>>>()),
                       bagel_determinant_map_(std::make_shared<std::map< std::string, std::shared_ptr<const Determinants>>>()) {};
    ~CI_Type_Converter(){};
   
    void add_civec( std::shared_ptr<const Civec> civec, int state_num );
   
  };
}
}
#endif
