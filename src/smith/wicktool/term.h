#ifndef __SRC_SMITH_WICKTOOL_BRAKET_INIT_H
#define __SRC_SMITH_WICKTOOL_BRAKET_INIT_H

#include <src/smith/wicktool/wickutils.h>

template<typename DataType>
class BraKet_Init {

     public :
       const std::vector<std::string> op_list;
       const DataType factor;
       const int bra_num_;
       const int ket_num_;
       const std::string type ; // should be "ci_deriv" or "full" 
       const std::string multiop;

       BraKet_Init( std::pair< std::vector<std::string>, DataType > BraKet_info, int bra_num, int ket_num, std::string type_in ) :
                  op_list(BraKet_info.first), factor(BraKet_info.second), bra_num_(bra_num), ket_num_(ket_num), type(type_in), multiop("asd") {} ;
                 // multiop(std::accumulate(op_list.begin(), op_list.end(), std::string(""))) {};
       ~BraKet_Init(){};
};
#endif
