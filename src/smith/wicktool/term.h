#ifndef __SRC_SMITH_WICKTOOL_TERM_H
#define __SRC_SMITH_WICKTOOL_TERM_H
#include <src/smith/wicktool/expression.h>
#include <src/smith/wicktool/braket.h>
#include <src/smith/tensor.h>
#include <src/smith/wicktool/states_info.h>

template<class DataType> 
class Term_Info {

     public :
       const std::vector<std::string> op_list;
       const DataType factor;
       const std::string Bra_name;
       const std::string Ket_name;
       const std::string type ; // should be "ci_deriv" or "full" 
       const std::string term_name;

     public :
       Term_Info( std::pair< std::vector<std::string>, DataType > BraKet_info, 
                  std::string Bra_name_in, std::string Ket_name_in, std::string type_in ) :
                  op_list(BraKet_info.first), factor(BraKet_info.second), Bra_name(Bra_name_in), Ket_name(Ket_name_in),
                  type(type_in) {};
       ~Term_Info(){};
};


//template<class DataType> 
//class Term {
//     using BraKet_Info = Term_Info<DataType> ;
//
//     public :
//       const std::vector<std::string> ;
//       const DataType factor;
//       const std::string Bra_name;
//       const std::string Ket_name;
//       const std::string type ; // should be "ci_deriv" or "full" 
//       const std::string term_name;

//     public :
//       Term_Info( std::pair< std::vector<std::string>, DataType > BraKet_info, 
//                  std::string Bra_name_in, std::string Ket_name_in, std::string type_in ) :
//                  op_list(BraKet_info.first), factor(BraKet_info.second), Bra_name(Bra_name_in), Ket_name(Ket_name_in),
//                  type(type_in) {};
//       ~Term_Info(){};
//};
#endif
