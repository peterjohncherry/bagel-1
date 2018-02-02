#ifndef __SRC_PROP_PROPTOOL_TASKTRANSLATOR_SYSTEMCOMPUTER_H
#define __SRC_PROP_PROPTOOL_TASKTRANSLATOR_SYSTEMCOMPUTER_H

#include <src/prop/proptool/proputils.h> 
#include <src/prop/proptool/algebraic_manipulator/system_info.h>
#include <src/prop/proptool/task_translator/equation_computer.h>

namespace  System_Computer { 
template<class DataType> 
class System_Computer {

  public:
    System_Computer();
   ~System_Computer(){};

    void build_equation_computer(std::string equation_name ); 
    void build_expression_computer( std::string expression_name ) ;
    void build_tensop( std::string tensop_name ) ;

};
}
#endif
