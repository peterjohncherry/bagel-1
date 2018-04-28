#ifndef __SRC_PROP_PROPTOOL_a_contrib_info_H
#define __SRC_PROP_PROPTOOL_a_contrib_info_H

#include <src/prop/proptool/proputils.h>

using namespace WickUtils;

class AContribInfo_Base {

  public :
    std::string name_;
    int total_uses_;
    int remaining_uses_;

    AContribInfo_Base( std::string name ) : name_(name), total_uses_(1), remaining_uses_(1) {};

    ~AContribInfo_Base(){};

    std::string name() {return name_ ;}

    int total_uses() const { return total_uses_ ; }
    int remaining_uses() const { return remaining_uses_ ; }

    void increase_total_uses() { total_uses_+=1; }
    void increase_remaining_uses() { remaining_uses_+=1; }
    void decrease_total_uses() { total_uses_-=1; }
    void decrease_remaining_uses() { remaining_uses_-=1; }


    virtual std::vector<int> id_order(int qq) { assert(false); std::vector<int> dummy; return dummy; };
    virtual std::vector<std::vector<int>> id_orders() { assert(false); std::vector<std::vector<int>> dummy ;return dummy; };

    virtual void add_id_order( std::vector<int>& new_id_order ) { assert(false); } 
    virtual std::pair<double,double> factor(int qq) { assert(false); std::pair<double,double> dummy; return dummy; };

    virtual void add_factor( std::pair<double,double>& new_factor ) {};
    virtual void combine_factors( int qq, std::pair<double,double>& new_factor ) {}; 


};

template<typename DataType>
class AContribInfo_Full : public AContribInfo_Base {

  public :
    std::vector<std::vector<int>> id_orders_;
    std::vector<std::pair<double,double>> factors;

    AContribInfo_Full( std::string name,  std::vector<int>& id_order , std::pair<double,double> factor ):
                       AContribInfo_Base(name),  factors(std::vector<std::pair<double,double>>(1,factor)),
                       id_orders_(std::vector<std::vector<int>>(1,id_order)) {};
    ~AContribInfo_Full(){};

    std::pair<double,double> factor(int qq) {return factors[qq]; };
    void add_factor( std::pair<double,double>& new_factor ) { factors.push_back(new_factor); };
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors[qq].first += new_factor.first;  
      factors[qq].second += new_factor.second;  
    }

    std::vector<int> id_order(int qq) { return id_orders_[qq]; };
    std::vector<std::vector<int>> id_orders() { return id_orders_; };

    void add_id_order( std::vector<int>& new_id_order ) {  id_orders_.push_back(new_id_order); }  
};

#endif
