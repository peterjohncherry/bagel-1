#ifndef __SRC_PROP_PROPTOOL_CONSTRAINTS_H
#define __SRC_PROP_PROPTOOL_CONSTRAINTS_H

#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include <map> 
#include <numeric> 
#include <algorithm>
#include <functional>
#include <utility> 

class Constraint { 
  public :

   std::string name_;

   Constraint( std::string name ) : name_(name) {}  
   ~Constraint(){}
   
   virtual bool apply_constraint( const std::vector<std::string>& rngs ) = 0 ; 

};

class Constraint_Dynamic : public Constraint { 

  public:
   std::function< bool (const std::vector<std::string>& )> constraint_function_;

   Constraint_Dynamic( std::string name, std::function< bool(const std::vector<std::string>& )>& constraint_function ) :
                       Constraint( name ), constraint_function_(constraint_function) {} 

   ~Constraint_Dynamic(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs ) { return constraint_function_(rngs); }   

}; 


class Constraint_NotAllAct : public Constraint { 

  public : 

   Constraint_NotAllAct( std::string name ) : Constraint(name) {}  
   ~Constraint_NotAllAct(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs );  

};
#endif
