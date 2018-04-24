#ifndef __SRC_PROP_PROPTOOL_a_contrib_info_H
#define __SRC_PROP_PROPTOOL_a_contrib_info_H

#include <src/prop/proptool/proputils.h>

using namespace WickUtils;

template<typename DataType>
class AContribInfo {

  public :
    std::string name_;
    int total_uses_;
    int remaining_uses_;

    AContribInfo( std::string name ) : name_(name), total_uses_(1), remaining_uses_(1) {};

    ~AContribInfo(){};

    std::string name() {return name_ ;}

    int total_uses() const { return total_uses_ ; }
    int remaining_uses() const { return remaining_uses_ ; }

    void increase_total_uses() { total_uses_+=1; }
    void increase_remaining_uses() { remaining_uses_+=1; }
    void decrease_total_uses() { total_uses_-=1; }
    void decrease_remaining_uses() { remaining_uses_-=1; }

    virtual std::pair<DataType,DataType> factor(int qq) = 0;
    virtual std::pair<DataType,DataType> factor(int qq, int rr) = 0;
    virtual void add_factor( std::pair<DataType,DataType>& new_factor ) = 0;
    virtual void add_factor( int qq, std::pair<DataType,DataType>& new_factor ) = 0;
    virtual void combine_factors( int qq, std::pair<DataType,DataType>& new_factor ) = 0;
    virtual void combine_factors( int qq, int rr, std::pair<DataType,DataType>& new_factor ) = 0;

    virtual std::vector<std::vector<int>>  id_orders() = 0; // TODO change this to ptr
    virtual std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>> pid_orders() = 0; 
    virtual std::vector<std::shared_ptr<std::vector<int>>> aid_orders() = 0;

    virtual std::vector<int> id_order(int qq) = 0; // TODO change this to ptr
    virtual std::shared_ptr<std::vector<int>> aid_order(int qq) = 0;
    virtual std::shared_ptr<std::vector<int>> pid_order(int qq, int rr ) = 0; 
    virtual std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> aid_pid_orders(int qq) = 0;

    virtual void add_id_order(  std::vector<int>& new_id_order) = 0; // TODO change this to ptr
    virtual void add_aid_order( std::vector<int>& new_aid_order ) = 0;
    virtual void add_pid_order( int qq,  std::vector<int>& new_pid_order ) = 0; 

};

template<typename DataType>
class AContribInfo_Full : public AContribInfo<DataType> {

  public :
    std::vector<std::vector<int>> id_orders_;
    std::vector<std::pair<DataType,DataType>> factors;

    AContribInfo_Full( std::string name,  std::vector<int>& id_order , std::pair<DataType,DataType> factor ):
                  AContribInfo<DataType>(name),  factors(std::vector<std::pair<DataType,DataType>>(1,factor)),
                  id_orders_(std::vector<std::vector<int>>(1,id_order)) {};
    ~AContribInfo_Full(){};

    std::pair<DataType,DataType> factor(int qq) {return factors[qq]; };
    void add_factor( std::pair<DataType,DataType>& new_factor ) { factors.push_back(new_factor); };
 
    void combine_factors( int qq, std::pair<DataType,DataType>& new_factor )  {
      factors[qq].first += new_factor.first;  
      factors[qq].second += new_factor.second;  
    }

    std::vector<int> id_order(int qq) { return id_orders_[qq]; };
    std::vector<std::vector<int>> id_orders() { return id_orders_; };

    std::pair<DataType,DataType> factor(int qq, int rr) {
      throw std::logic_error( "should not call from AContribInfo_Full X0" ); return std::make_pair(-1.0,-1.0);
    }

    void add_factor( int qq, std::pair<DataType,DataType>& new_factor )  {
      throw std::logic_error( "should not call from AContribInfo_Full X1" );
    }

    void combine_factors( int qq, int rr, std::pair<DataType,DataType>& new_factor )  {
      throw std::logic_error( "should not call from AContribInfo_Full X2" );
    }

    std::shared_ptr<std::vector<int>> aid_order(int qq) { 
        throw std::logic_error( " should not call aid_order from AContribInfo_Full ") ; 
      return std::make_shared<std::vector<int>>();
    };

    std::shared_ptr<std::vector<int>> pid_order( int qq, int rr ) {
        throw std::logic_error( " should not call pid_order from AContribInfo_Full ") ; 
      return std::make_shared<std::vector<int>>();
    };

    std::vector<std::shared_ptr<std::vector<int>>> aid_orders() {
        throw std::logic_error( "should not call aid_orders from AContribInfo_Full ") ; 
        std::vector<std::shared_ptr<std::vector<int>>> dummy;
      return dummy;
    }

    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> aid_pid_orders(int qq) {
       throw std::logic_error( "should not call aid_pid_orders from AContribInfo_Full ") ; 
       std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> dummy; 
      return dummy;
    }

    std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>> pid_orders() {
       throw std::logic_error( "should not call pid_orders from AContribInfo_Full ") ; 
       std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>> dummy; 
      return dummy;
    }

    void add_id_order( std::vector<int>& new_id_order ) {  id_orders_.push_back(new_id_order); }  
    void add_aid_order( std::vector<int>& new_aid_order ) { throw std::logic_error( "should not call add_aid_order from AContribInfo_Full ") ; }  
    void add_pid_order( int qq, std::vector<int>& new_pid_order ) { throw std::logic_error( "should not call add_pid_order from AContribInfo_Full ") ; }  
};

template<typename DataType>
class AContribInfo_ExcDeriv : public AContribInfo<DataType> {

  public :
    std::vector<std::shared_ptr<std::vector<int>>> aid_orders_;
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>> pid_orders_;
    std::vector<std::vector<std::pair<DataType,DataType>>> pid_factors_;

    AContribInfo_ExcDeriv( std::string name,  std::vector<int>& aid_order, std::vector<int>& pid_order,
                  std::pair<DataType,DataType> factor ): AContribInfo<DataType>(name),
                  aid_orders_(std::vector<std::shared_ptr<std::vector<int>>>(1,std::make_shared<std::vector<int>>(aid_order))),
                  pid_orders_(std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>>(1,
                               std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>> (1, std::make_shared<std::vector<int>>(pid_order)))),
                  pid_factors_(std::vector<std::vector<std::pair<DataType,DataType>>>( 1, std::vector<std::pair<DataType,DataType>>(1, factor))) 
                  {};

    ~AContribInfo_ExcDeriv(){};

    std::pair<DataType,DataType> factor(int qq) {
      throw std::logic_error( "should not call from AContribInfo_ExcDeriv X0" ); return std::make_pair(-1.0,-1.0);
    }

    void add_factor( std::pair<DataType,DataType>& new_factor )  {
      throw std::logic_error( "should not call from AContribInfo_ExcDeriv X1" );
    }

    void combine_factors( int qq, std::pair<DataType,DataType>& new_factor )  {
      throw std::logic_error( "should not call from AContribInfo_ExcDeriv X2" );
    }

    std::vector<int> id_order(int qq) {
        throw std::logic_error( " should not call id_order from AContribInfo_Exc ") ; 
      return *(aid_orders_[qq]);
    };
 
    std::vector<std::vector<int>> id_orders() {
        throw std::logic_error( " should not call id_orders from AContribInfo_Exc ") ; 
       std::vector<std::vector<int>> dummy(0); 
      return dummy;
    };
 
    std::pair<DataType,DataType> factor(int qq, int rr) {
       return pid_factors_[qq][rr];
    }

    void add_factor( int qq, std::pair<DataType,DataType>& new_factor )  {
      if ( pid_factors_.size() < qq ) {
        pid_factors_.push_back( std::vector<std::pair<DataType,DataType>>(1, new_factor) );
      } else {  
        pid_factors_[qq].push_back(new_factor);  
      } 
    }

    void combine_factors( int qq, int rr, std::pair<DataType,DataType>& new_factor )  {
        pid_factors_[qq][rr].first +=new_factor.first;  
        pid_factors_[qq][rr].second +=new_factor.second;  
    }

    std::shared_ptr<std::vector<int>> aid_order(int qq) { return aid_orders_[qq]; }
    std::shared_ptr<std::vector<int>> pid_order(int qq, int rr) { return pid_orders_[qq]->at(rr); }
   
    std::vector<std::shared_ptr<std::vector<int>>> aid_orders() { return aid_orders_; }
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>>> pid_orders() { return pid_orders_; }

    void add_id_order( std::vector<int>& new_id_order ) { throw std::logic_error( "should not call add_id_order from AContribInfo_Exc ") ; }  

    void add_aid_order( std::vector<int>& new_aid_order ) { aid_orders_.push_back(std::make_shared<std::vector<int>>(new_aid_order)); }  

    void add_pid_order( int qq,  std::vector<int>& new_pid_order ) {
      if ( pid_orders_.size() < qq ) {
        pid_orders_.push_back(std::make_shared<std::vector<std::shared_ptr<std::vector<int>>>>( 1, std::make_shared<std::vector<int>>(  new_pid_order)));
      } else {
        pid_orders_[qq]->push_back(std::make_shared<std::vector<int>>(new_pid_order));
      }
    }
 
    std::shared_ptr<std::vector<std::shared_ptr<std::vector<int>>>> aid_pid_orders(int qq) {
      return pid_orders_[qq]; 
    }
  
};

#endif
