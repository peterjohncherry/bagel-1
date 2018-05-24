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

    virtual std::vector<int>& id_order(int qq) = 0;
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
    std::shared_ptr<std::vector<int>> post_contraction_reordering() { throw std::logic_error("Do not call post_contraction_reordering() from AContribInfo_Base"); return std::make_shared<std::vector<int>>(0); }
    
    virtual
    std::vector<int> gamma_contraction_pos() {  throw std::logic_error( "Should not access gamma_contraction_pos from AContribInfo_Base" );  return std::vector<int>(0); }

    virtual 
    std::shared_ptr<std::map<std::vector<int>, std::shared_ptr<std::map<std::string , std::shared_ptr<AContribInfo_Base>>>>> gamma_pos_map(){ 
      throw std::logic_error( "Should not access gamma_contraction_pos from AContribInfo_Base" );
      return std::make_shared<std::map<std::vector<int>, std::shared_ptr<std::map<std::string , std::shared_ptr<AContribInfo_Base>>>>>();
    }

    virtual std::shared_ptr<std::vector<std::string>> post_gamma_contraction_ranges() { 
       throw std::logic_error( "Should not access post_gamma_contraction_ranges from AContribInfo_Base" );
    return std::make_shared<std::vector<std::string>>(0); }

    virtual  std::shared_ptr<std::vector<std::string>> a_block_ranges() {
       throw std::logic_error( "Should not access a_block_ranges from AContribInfo_Base" );
    return std::make_shared<std::vector<std::string>>(0); }

    virtual 
    std::shared_ptr<std::vector<int>> target_block_positions() { 
      throw std::logic_error( "Should not access target_block_positions from AContribInfo_Base" );
    return std::make_shared<std::vector<int>>(); }
};

template<typename DataType>
class AContribInfo_Full : public AContribInfo_Base {

  public :
    std::vector<std::vector<int>> id_orders_;
    std::vector<std::pair<double,double>> factors;
    std::shared_ptr<std::vector<std::string>> post_reorder_rngs_;

    AContribInfo_Full( std::string name,  std::vector<int>& id_order , std::pair<double,double> factor ):
                       AContribInfo_Base(name),  factors(std::vector<std::pair<double,double>>(1,factor)),
                       id_orders_(std::vector<std::vector<int>>(1,id_order)) {};

    AContribInfo_Full( std::string name,  std::vector<int>& id_order, std::shared_ptr<std::vector<std::string>> post_reorder_rngs, std::pair<double,double> factor ):
                       AContribInfo_Base(name),  factors(std::vector<std::pair<double,double>>(1,factor)),
                       id_orders_(std::vector<std::vector<int>>(1,id_order)), post_reorder_rngs_(post_reorder_rngs) {};

    ~AContribInfo_Full(){};

    std::pair<double,double> factor(int qq) {return factors[qq]; };

    void add_factor( std::pair<double,double>& new_factor ) { factors.push_back(new_factor); };
    void add_id_order( std::vector<int>& new_id_order ) { id_orders_.push_back(new_id_order);}
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors[qq].first += new_factor.first;  
      factors[qq].second += new_factor.second;  
    }

    std::vector<int>& id_order(int qq) { return id_orders_[qq]; };
    std::vector<std::vector<int>> id_orders() { return id_orders_; };
 
    std::shared_ptr<std::vector<std::string>> a_block_ranges() { return post_reorder_rngs_ ; } 

};

template<typename DataType>
class AContribInfo_OrbExcDeriv : public AContribInfo_Base {

  public :
    std::string target_block_name_;
    std::vector<std::pair<double,double>> factors_;
    std::shared_ptr<std::vector<int>> target_block_positions_;
    std::shared_ptr<std::vector<std::string>> post_gamma_contraction_ranges_;
    std::shared_ptr<std::vector<int>> post_contraction_reordering_;

    // key : rearrangement after gamma contraction
    // key : rearrangement before gamma contraction
    std::shared_ptr<std::map<std::vector<int>, std::shared_ptr<std::map<std::string , std::shared_ptr<AContribInfo_Base>>>>>  gamma_pos_map_;


    AContribInfo_OrbExcDeriv( std::string ablock_name, std::string target_block_name, 
                              std::shared_ptr<std::vector<int>> target_block_positions, std::shared_ptr<std::vector<std::string>> post_gamma_contraction_ranges ):
                              AContribInfo_Base(ablock_name), target_block_name_( target_block_name ), target_block_positions_(target_block_positions),
                              post_gamma_contraction_ranges_(post_gamma_contraction_ranges) { 
                                gamma_pos_map_=std::make_shared<std::map<std::vector<int>, std::shared_ptr<std::map<std::string,std::shared_ptr<AContribInfo_Base>>>>>();
                                post_contraction_reordering_ = get_ascending_order( *target_block_positions_ );
                                assert(post_gamma_contraction_ranges_->size() == post_contraction_reordering_->size());
                              };

    ~AContribInfo_OrbExcDeriv(){};

    std::pair<double,double> factor( int qq ) { return factors_[qq]; } ;

    void add_reordering( std::string ablock_name,  std::vector<int>& gamma_contraction_pos,
                         std::vector<int>& pre_contraction_reordering,
                         std::shared_ptr<std::vector<std::string>> pre_contraction_ranges, std::pair<double,double> factor ) {

      auto ablock_map_loc = gamma_pos_map_->find( gamma_contraction_pos ); 
      if ( ablock_map_loc == gamma_pos_map_->end() ) {

        auto ablock_map = std::make_shared<std::map< std::string , std::shared_ptr<AContribInfo_Base> >>(); 
        auto a_info = std::make_shared<AContribInfo_Full<DataType>>( ablock_name, pre_contraction_reordering, pre_contraction_ranges, factor );
        ablock_map->emplace( ablock_name, a_info);
        gamma_pos_map_->emplace( gamma_contraction_pos, ablock_map );
    
      } else {

        std::shared_ptr<std::map< std::string, std::shared_ptr<AContribInfo_Base>>> ablock_map = ablock_map_loc->second;
        auto ablock_loc = ablock_map->find(ablock_name);
        if ( ablock_loc == ablock_map->end() ) {

          auto a_info = std::make_shared<AContribInfo_Full<DataType>>( ablock_name, pre_contraction_reordering, pre_contraction_ranges, factor );
          ablock_map->emplace( ablock_name, a_info );
          gamma_pos_map_->emplace( gamma_contraction_pos, ablock_map );

        } else { 
          std::shared_ptr<AContribInfo_Full<DataType>> AInfo = std::dynamic_pointer_cast<AContribInfo_Full<DataType>>(ablock_loc->second);
 
          for ( int qq = 0 ; qq != AInfo->id_orders().size(); qq++ ) {
            if( pre_contraction_reordering == AInfo->id_order(qq) ){
              AInfo->combine_factors( qq, factor );
              AInfo->remaining_uses_ += 1;
              AInfo->total_uses_ += 1;
              break;
          
            } else if ( qq == AInfo->id_orders().size()-1) {
              AInfo->add_id_order(pre_contraction_reordering);
              AInfo->add_factor(factor);
            }
          }
        }
      }
    }
 
    void combine_factors( int qq, std::pair<double,double>& new_factor )  {
      factors_[qq].first += new_factor.first;  
      factors_[qq].second += new_factor.second;  
    }

    std::string target_block_name() { return target_block_name_; } 
     
    std::shared_ptr<std::map<std::vector<int>, std::shared_ptr<std::map<std::string , std::shared_ptr<AContribInfo_Base>>>>> gamma_pos_map(){ return gamma_pos_map_ ;}
 
    std::shared_ptr<std::vector<std::string>> post_gamma_contraction_ranges() { return post_gamma_contraction_ranges_; }
    std::shared_ptr<std::vector<int>> post_contraction_reordering() { return post_contraction_reordering_; }
    std::shared_ptr<std::vector<int>> target_block_positions() { return target_block_positions_; }
    std::vector<int>& id_order(int qq) { return *target_block_positions_; };

};
#endif
