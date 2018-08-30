#include <bagel_config.h>
#include <src/prop/proptool/tensor_and_ci_lib/ci_type_converter.h>
#include <src/prop/proptool/proputils.h>

using namespace std;
using namespace bagel;
using namespace WickUtils;

//#define __DEBUG_PROPTOOL_CI_TYPE_CONVERTER
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Need different one for rel and non-rel
template<typename DataType>
void CI_Type_Converter::CI_Type_Converter<DataType>::add_civec(std::shared_ptr<const Civec> civec, int state_num )  {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __DEBUG_PROPTOOL_CI_TYPE_CONVERTER
cout << "CI_Type_Converter::CI_Type_Converter" << endl;
#endif //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //assert( range_conversion_map_ != NULL ); 

  std::string civec_name = get_civec_name( state_num, civec->det()->norb(), civec->det()->nelea(), civec->det()->neleb() );

  if ( civec_data_map_->find(civec_name) == civec_data_map_->end()){
 
    vector<IndexRange> civec_idxrng(1, *(range_conversion_map_->at(civec_name)) );

    shared_ptr<Tensor> civec_tensor = make_shared<Tensor>( civec_idxrng );
    civec_tensor->allocate();
    civec_tensor->zero();
    size_t coeff_pos_offset = 0;
    for ( Index idx_block : civec_idxrng[0].range() ){
       unique_ptr<DataType[]> civec_block(new DataType[idx_block.size()]);
       std::fill_n(civec_block.get(), idx_block.size(), DataType(0.0) );
       copy_n( civec->data() + coeff_pos_offset, idx_block.size(), civec_block.get());
       civec_tensor->put_block(civec_block, vector<Index>({ idx_block })) ;
       coeff_pos_offset += idx_block.size();
    }
    civec_data_map_->emplace( civec_name, civec_tensor);
  } else { 
    throw logic_error("Already have civec " +  civec_name + " in the map! Aborting!!" );
  } 
  return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template class CI_Type_Converter::CI_Type_Converter<double>;
template class CI_Type_Converter::CI_Type_Converter<std::complex<double>>;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
