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


    virtual std::vector<int> id_order(int qq) { assert(false); return std::vector<int>(0); }
    virtual std::vector<std::vector<int>> id_orders() { assert(false); return std::vector<std::vector<int>>(0); }

    virtual void add_id_order( std::vector<int>& new_id_order ) { assert(false); } 
    virtual std::pair<double,double> factor(int qq) { assert(false); return std::pair<double,double>(0.0,0.0); }

    virtual void add_factor( std::pair<double,double>& new_factor ) { throw std::logic_error("do not call add_factor from AContribInfo_Base"); }
    virtual void combine_factors( int qq, std::pair<double,double>& new_factor ) {throw std::logic_error("do not call combine_factors from AContribInfo_Base");  }

    virtual
    void add_pre_contraction_reordering( std::vector<int>& new_id_order ) {
      throw std::logic_error("do not call add_pre_contraction_reordering() from AContribInfo_Base"); }

    virtual std::string target_block_name() { throw  std::logic_error("do not call target_block_name() from AContribInfo_Base"); return "ERROR!"; } 

    virtual
    std::vector<std::vector<int>> pre_contraction_reorderings() {
      throw  std::logic_error("do not call pre_contraction_reorderings() from AContribInfo_Base");return std::vector<std::vector<int>>(0); }

    virtual
    std::vector<int> post_contraction_reordering() { throw std::logic_error("Do not call post_contraction_reordering() from AContribInfo_Base");return std::vector<int>(0); }
    
    virtual
    std::vector<int> gamma_contraction_pos() {  throw std::logic_error( "Should not access gamma_contraction_pos from AContribInfo_Base" );  return std::vector<int>(0); }

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
    void add_id_order( std::vector<int>& new_id_order ) { id_orders_.push_back(new_id_order);}
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors[qq].first += new_factor.first;  
      factors[qq].second += new_factor.second;  
    }

    std::vector<int> id_order(int qq) { return id_orders_[qq]; };
    std::vector<std::vector<int>> id_orders() { return id_orders_; };

};


template<typename DataType>
class AContribInfo_Part : public AContribInfo_Base {

  public :
    std::vector<std::vector<int>> id_orders_;
    std::vector<std::pair<double,double>> factors;

    AContribInfo_Part( std::string name, std::vector<int>& id_order, std::pair<double,double> factor ):
                       AContribInfo_Base(name),  factors(std::vector<std::pair<double,double>>(1,factor)), 
                       id_orders_(std::vector<std::vector<int>>(1,id_order)) {};
    ~AContribInfo_Part(){};

    std::pair<double,double> factor(int qq) {return factors[qq]; };

    void add_factor( std::pair<double,double>& new_factor ) { factors.push_back(new_factor); };
    void add_id_order( std::vector<int>& new_id_order ) { id_orders_.push_back(new_id_order);}
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors[qq].first += new_factor.first;  
      factors[qq].second += new_factor.second;  
    }

    std::vector<int> id_order(int qq) { return id_orders_[qq]; };
    std::vector<std::vector<int>> id_orders() { return id_orders_; };

};


template<typename DataType>
class AContribInfo_OrbExcDeriv : public AContribInfo_Base {

  public :
    std::string target_block_name_;
    std::vector<int> gamma_contraction_pos_;
    std::vector<std::pair<double,double>> factors_;

    // key : rearrangement after gamma contraction
    // key : rearrangement before gamma contraction
    std::map< std::vector<int> , std::shared_ptr<AContribInfo_Part<DataType>> > reorderings_map_;


    AContribInfo_OrbExcDeriv( std::string name, std::string target_block_name, std::vector<int>& gamma_contraction_pos,
                              std::vector<int>& pre_contraction_reordering, std::vector<int>& post_contraction_reordering,
                              std::pair<double,double> factor ):
                              AContribInfo_Base(name), target_block_name_( target_block_name ),
                              factors_(std::vector<std::pair<double,double>>(1,factor)) {
                                reorderings_map_.emplace( post_contraction_reordering, std::make_shared<AContribInfo_Part<DataType>>(name, pre_contraction_reordering, factor) );
                              };
    ~AContribInfo_OrbExcDeriv(){};

    std::pair<double,double> factor(int qq) {return factors_[qq]; };

    void add_reordering( std::vector<int>& post_contraction_reordering,  std::vector<int>& pre_contraction_reordering,
                         std::pair<double,double> factor ) {
      auto reorderings_map_loc = reorderings_map_.find( post_contraction_reordering );
      if ( reorderings_map_loc == reorderings_map_.end() ) {
        reorderings_map_.emplace( post_contraction_reordering, std::make_shared<AContribInfo_Part<DataType>>(name_, pre_contraction_reordering, factor));
      } else {
        std::shared_ptr<AContribInfo_Part<DataType>> AInfo = reorderings_map_loc->second;
        for ( int qq = 0 ; qq != AInfo->id_orders().size(); qq++ ) {
          if( pre_contraction_reordering == AInfo->id_orders_[qq] ){
            AInfo->combine_factors( qq, factor );
            AInfo->remaining_uses_ += 1;
            AInfo->total_uses_ += 1;
            break;
     
          } else if ( qq == AInfo->id_orders().size()-1) {
            AInfo->id_orders().push_back(pre_contraction_reordering);
            AInfo->add_factor(factor);
          }
        }
      }
    }
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors_[qq].first += new_factor.first;  
      factors_[qq].second += new_factor.second;  
    }

    std::string target_block_name() { return target_block_name_; } 
    std::vector<int> gamma_contraction_pos() { return gamma_contraction_pos_; } // defined w.r.t. order after pre_contraction_reordering

};


#endif
