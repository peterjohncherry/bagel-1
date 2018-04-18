#ifndef __SRC_PROP_PROPTOOL_SYMMOPS_H
#define __SRC_PROP_PROPTOOL_SYMMOPS_H

#include <memory>
#include <vector>
#include <tuple>
#include <string>
#include <map> 
#include <numeric> 
#include <algorithm>
#include <functional>
#include <utility> 

class Transformation {
  public :
 
    const std::string name_; 
    
    Transformation( const std::string name ): name_(name) {}
    ~Transformation(){}

    virtual std::vector<std::string> transform( const std::vector<std::string>& ranges ) = 0 ;

    virtual std::pair<double,double> factor( const std::vector<std::string>& ranges ) = 0 ;

    virtual std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& ranges ) = 0 ;

};
  

//For testing purposes; does nothing
class Transformation_ID : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, 1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_ID( const std::string name ): Transformation( name ),
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 0, 1, 2, 3} ) ), 
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 0, 1, 2, 3} ) ) {} 

    ~Transformation_ID(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ) { std::vector<std::string> rngs_new = rngs; return rngs_new; } ; 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ) { std::vector<std::string> rngs_new = rngs; return rngs_new; } ; 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};
 
class Transformation_UserDefined : public  Transformation {
  public : 

    std::vector<int> idxs_trans_;
    std::pair<double,double> factor_;
    std::function<std::vector<std::string>(const std::vector<std::string>& )> transform_;

    Transformation_UserDefined( const std::string name, const std::vector<int>& idxs_trans,  const std::pair<double,double>& factor, 
                                std::function<std::vector<std::string>(const std::vector<std::string>& )>& transform ):
                                Transformation( name ), idxs_trans_(idxs_trans), factor_(factor), transform_(transform) {}
    ~Transformation_UserDefined(){}

    std::vector<std::string> transform( const std::vector<std::string>& ranges ) { return transform_(ranges); };

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_ ; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& ranges ) { return std::make_shared<std::vector<int>>(idxs_trans_); }  ;

};
 
class Transformation_Hermitian : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, -1.0};

    Transformation_Hermitian( const std::string name ): Transformation( name )  {}
    ~Transformation_Hermitian(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) {
                                                  std::vector<int> new_order(rngs.size()) ;
                                                  std::iota( new_order.begin(), new_order.end(), 0 );
                                                  std::reverse( new_order.begin(), new_order.end());                     
                                                  return std::make_shared<std::vector<int>>( new_order ) ;
                                                } 

};
  
class Transformation_Spinflip : public Transformation {
  public :

    const std::pair<double,double> factor_ = {-1.0, 1.0};

    Transformation_Spinflip( const std::string name ): Transformation( name )  {}
    ~Transformation_Spinflip(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) {
                                                  std::vector<int> new_order(rngs.size()) ;
                                                  std::iota( new_order.begin(), new_order.end(), 0 );
                                                  std::reverse( new_order.begin(), new_order.end());                     
                                                  return std::make_shared<std::vector<int>>( new_order ) ;
                                                } 

};
 
class Transformation_1032 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, 1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_1032( const std::string name ): Transformation( name ),
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 1, 0, 3, 2} ) ), 
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 1, 0, 3, 2} ) ) {} 

    ~Transformation_1032(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};
       
class Transformation_2301 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, 1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_2301( const std::string name ): Transformation( name ),
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 2, 3, 0, 1 } ) ),
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 2, 3, 0, 1 } ) ) {} 
    ~Transformation_2301(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};
 
class Transformation_2103 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, 1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_2103( const std::string name ): Transformation( name ),
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 2, 1, 0, 3 } ) ),
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 2, 1, 0, 3 } ) ) {} 

    ~Transformation_2103(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};
 

class Transformation_3012 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, -1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_3012( const std::string name ): Transformation( name ) ,
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 3, 0, 1, 2 } ) ) ,
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 1, 2, 3, 0 } ) ) {} 

    ~Transformation_3012(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};
 

class Transformation_0321 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, -1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_0321( const std::string name ): Transformation( name ) ,
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 0, 3, 2, 1} ) ),
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 0, 3, 2, 1} ) ) {} 

    ~Transformation_0321(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};

class Transformation_1230 : public Transformation {
  public :

    const std::pair<double,double> factor_ = {1.0, -1.0};
    std::shared_ptr<std::vector<int>> idxs_trans_; 
    std::shared_ptr<std::vector<int>> idxs_trans_inverse_; 
 
    Transformation_1230( const std::string name ): Transformation( name ) ,
           idxs_trans_( std::make_shared<std::vector<int>>( std::vector<int> { 1, 2, 3, 0 } ) ) ,
           idxs_trans_inverse_( std::make_shared<std::vector<int>>( std::vector<int> { 3, 0, 1, 2 } ) ) {} 

    ~Transformation_1230(){}

    std::vector<std::string> transform( const std::vector<std::string>& rngs ); 
    std::vector<std::string> inverse_transform( const std::vector<std::string>& rngs ); 

    std::pair<double,double> factor( const std::vector<std::string>& ranges ) { return factor_; }

    std::shared_ptr<std::vector<int>> idxs_trans( const std::vector<std::string>& rngs ) { return idxs_trans_;} 

};

#endif
