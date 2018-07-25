#ifndef __SRC_PROP_PROPTOOL_CONSTRAINTS_H
#define __SRC_PROP_PROPTOOL_CONSTRAINTS_H

#include <memory>
#include <vector>
#include <string>
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
   
   Constraint_NotAllAct() : Constraint("Not_All_Act") {}  
   ~Constraint_NotAllAct(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs );  

};

class Constraint_Spin_Neutral_Normal_Order : public Constraint { 

  public : 

   Constraint_Spin_Neutral_Normal_Order() : Constraint("Spin_Neutral_Normal_Order") {}  
   ~Constraint_Spin_Neutral_Normal_Order(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs );  

};

class Constraint_All_Same_Spin : public Constraint { 

  public : 

   Constraint_All_Same_Spin() : Constraint("All_Same_Spin") {}  
   ~Constraint_All_Same_Spin(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs );  

};

class Constraint_Allowed_Blocks : public Constraint { 

  public : 
   
   std::vector<std::vector<std::string>> allowed_blocks_;

   Constraint_Allowed_Blocks() : Constraint("Allowed_Block_List") {};  

   Constraint_Allowed_Blocks( const std::vector<std::vector<std::string>>& allowed_blocks ) :
                              Constraint("Allowed_Block_List"), allowed_blocks_(allowed_blocks){};  

   ~Constraint_Allowed_Blocks(){}
   
   bool apply_constraint( const std::vector<std::string>& rngs );  

   void add_block( const std::vector<std::string>& block ){ allowed_blocks_.push_back(block); }

};



#endif
