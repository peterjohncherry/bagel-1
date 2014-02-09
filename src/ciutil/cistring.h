//
// BAGEL - Parallel electron correlation program.
// Filename: cistring.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef BAGEL_CIUTIL_STRINGSPACE_H
#define BAGEL_CIUTIL_STRINGSPACE_H

#include <vector>
#include <array>
#include <memory>
#include <bitset>
#include <algorithm>

#include <src/util/constants.h>
#include <src/parallel/staticdist.h>
#include <src/parallel/mpi_interface.h>
#include <src/util/serialization.h>
#include <src/ciutil/cistringmap.h>

namespace bagel {

// Contains all the strings and information for lexical ordering for one particular graph (set of strings)
//   comprised of three subgraphs (one each for RASI, RASII, RASIII)
class CIGraph {
  protected:
    size_t nele_;
    size_t norb_;
    size_t size_;
    std::vector<size_t> weights_;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & nele_ & norb_ & size_& weights_;
    }

  public:
    CIGraph() { }
    CIGraph(const size_t nele, const size_t norb);

    size_t& weight(const size_t i, const size_t j) { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }
    const size_t& weight(const size_t i, const size_t j) const { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }

    size_t size() const { return size_; }

    size_t lexical(const int& start, const int& fence, const std::bitset<nbit__>& abit) const {
      size_t out = 0;

      int k = 0;
      for (int i = start; i < fence; ++i)
        if (abit[i]) out += weight(i-start,k++);
      return out;
    }
};


class CIString_base {
  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) { }
  public:
    CIString_base() { }
    virtual ~CIString_base() { }
};


template<int N>
class CIString_base_impl : public CIString_base {
  protected:
    int norb_;
    int nele_;
    size_t offset_;
    std::vector<std::bitset<nbit__>> strings_;

    std::array<std::pair<int, int>, N> subspace_;
    std::array<std::shared_ptr<CIGraph>, N> graphs_;
    std::shared_ptr<const StaticDist> dist_;

    std::shared_ptr<StringMap> phi_;
    std::shared_ptr<StringMap> uncompressed_phi_;

    void init() {
      compute_strings();
      construct_phi();
    }
    virtual void compute_strings() = 0;
    virtual void construct_phi() = 0;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template <class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<CIString_base>(*this) << norb_ << nele_ << offset_ << strings_ << subspace_ << graphs_;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<CIString_base>(*this) >> norb_ >> nele_ >> offset_ >> strings_ >> subspace_ >> graphs_;
      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
    }

  public:
    CIString_base_impl() { }
    CIString_base_impl(std::initializer_list<size_t> args) {
      assert(args.size() == 2*N+1);
      auto iter = args.begin();
      for (int i = 0; i != N; ++i) {
        const int a = *iter++;
        const int b = *iter++;
        subspace_[i] = std::make_pair(a, b);
        graphs_[i] = std::make_shared<CIGraph>(a, b);
      }
      // setting to CIString_base
      offset_ = *iter++;
      norb_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.second; });
      nele_ = std::accumulate(subspace_.begin(), subspace_.end(), 0, [](int n, const std::pair<int, int>& i) { return n+i.first; });
      assert(iter == args.end());

      const size_t size = std::accumulate(graphs_.begin(), graphs_.end(), 1, [](size_t n, const std::shared_ptr<CIGraph>& i) { return n*i->size(); });
      dist_ = std::make_shared<StaticDist>(size, mpi__->size());
    }

    virtual ~CIString_base_impl() { }

    int nele() const { return nele_; }
    int norb() const { return norb_; }
    size_t size() const { return strings_.size(); }
    size_t offset() const { return offset_; }

    size_t key() const {
      size_t out = 0;
      for (int i = 0; i != N; ++i) {
        assert(subspace_[i].first < (1 << 6));
        out = (out << 6) + subspace_[i].first;
      }
      return out;
    }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__>& strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

    size_t size(const int& i) const { return graphs_[i]->size(); }

    std::shared_ptr<const StaticDist> dist() const { return dist_; }

    std::shared_ptr<const StringMap> phi() const { return phi_; }
    std::shared_ptr<const StringMap> uncompressed_phi() const { return uncompressed_phi_; }

    virtual size_t lexical_zero(const std::bitset<nbit__>& bit) const = 0;
    virtual size_t lexical_offset(const std::bitset<nbit__>& bit) const = 0;

    virtual bool contains(const std::bitset<nbit__>& bit) const = 0;
};


class RASString : public CIString_base_impl<3> {
  protected:
    void compute_strings() override;
    void construct_phi() override { /* TODO to be implemented */ }

    // helper functions
    int nholes(const std::bitset<nbit__>& bit) const {
      return subspace_[0].second - (bit & std::bitset<nbit__>((1ull << subspace_[0].second) - 1ull)).count();
    }
    int nparticles(const std::bitset<nbit__>& bit) const {
      return (bit & std::bitset<nbit__>(((1ull << subspace_[2].second) - 1ull) << (subspace_[0].second + subspace_[1].second))).count();
    }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<CIString_base_impl<3>>(*this);
    }

  public:
    RASString() { }
    RASString(const size_t nele1, const size_t norb1, const size_t nele2, const size_t norb2, const size_t nele3, const size_t norb3, const size_t offset = 0);

    int nholes() const { return subspace_[0].second - subspace_[0].first; }
    int nele2() const { return nele_ - subspace_[0].first - subspace_[2].first; }
    int nparticles() const { return subspace_[2].first; }

    size_t lexical_offset(const std::bitset<nbit__>& bit) const override {
      return lexical_zero(bit) + offset_;
    }

    size_t lexical_zero(const std::bitset<nbit__>& bit) const override {
      const size_t r1 = subspace_[0].second;
      const size_t r2 = subspace_[1].second;
      const size_t r3 = subspace_[2].second;

      const size_t n2 = graphs_[1]->size();
      const size_t n1 = graphs_[0]->size();

      size_t out = 0;
      out += graphs_[1]->lexical(r1, r1+r2, bit);
      out += n2 * graphs_[0]->lexical(0, r1, bit);
      out += n2 * n1 * graphs_[2]->lexical(r1+r2, r1+r2+r3, bit);

      return out;
    }

    template <int S> const std::pair<const int, const int> ras() const {
      static_assert(S == 0 || S == 1 || S == 2, "illegal call of RAString::ras");
      return std::get<S>(subspace_);
    }

    bool contains(const std::bitset<nbit__>& bit) const {
      assert(bit.count() == nele_);
      return nholes(bit) == nholes() && nparticles(bit) == nparticles();
    }
};


class FCIString : public CIString_base_impl<1> {
  protected:
    void compute_strings() override;
    void construct_phi() override;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<CIString_base_impl<1>>(*this);
    }

  public:
    FCIString() { }
    FCIString(const size_t nele1, const size_t norb1, const size_t offset = 0);

    size_t lexical(const std::bitset<nbit__>& bit) const {
      assert(contains(bit));
      return graphs_[0]->lexical(0, norb_, bit);
    }
    size_t lexical_offset(const std::bitset<nbit__>& bit) const override { return lexical(bit)+offset_; }
    size_t lexical_zero(const std::bitset<nbit__>& bit) const override { return lexical(bit); }

    bool contains(const std::bitset<nbit__>& bit) const { assert(bit.count() == nele_); return true; }


};

}

extern template class bagel::CIString_base_impl<1>;
extern template class bagel::CIString_base_impl<3>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::CIGraph)
BOOST_CLASS_EXPORT_KEY(bagel::CIString_base)
BOOST_CLASS_EXPORT_KEY(bagel::CIString_base_impl<1>)
BOOST_CLASS_EXPORT_KEY(bagel::CIString_base_impl<3>)
BOOST_CLASS_EXPORT_KEY(bagel::RASString)
BOOST_CLASS_EXPORT_KEY(bagel::FCIString)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<CIString_base, T>::value>::type> {
    typedef CIString_base type;
  };
}

#endif