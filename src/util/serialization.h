//
// BAGEL - Parallel electron correlation program.
// Filename: serialization.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
//

#ifndef __SRC_UTIL_SERIALIZATION_H
#define __SRC_UTIL_SERIALIZATION_H

#include <array>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <type_traits>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/archive/shared_ptr_helper.hpp>
#include <boost/property_tree/ptree_serialization.hpp>

// default implementation
namespace bagel {
  template <class T, typename = void>
  struct base_of {
    // always return non-const type
    typedef typename std::remove_cv<T>::type type;
  };
}

namespace boost {
  namespace serialization {

    template<class Archive, class T>
    inline void serialize(Archive& ar, std::shared_ptr<T>& t, const unsigned int version) {
      BOOST_STATIC_ASSERT(tracking_level<T>::value != track_never);
      split_free(ar, t, version);
    }

    template<class Archive, class T>
    inline void save(Archive& ar, const std::shared_ptr<T>& t, const unsigned int) {
      BOOST_STATIC_ASSERT((tracking_level<T>::value != track_never));
      T* ptr = t.get();
      ar << ptr;
    }

    // base, non-const classes
    template<class Archive, class T, typename = void>
    struct load_impl {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        BOOST_STATIC_ASSERT((tracking_level<T>::value != track_never));
        T* ptr;
        ar >> ptr;

        static std::unordered_map<T*, std::weak_ptr<T>> hash;
        if (hash[ptr].expired()) {
          t = std::shared_ptr<T>(ptr);
          hash[ptr] = t;
        } else {
          t = hash[ptr].lock();
        }
      }
    };

    // non-base non-const classes
    template<class Archive, class T>
    struct load_impl<Archive, T,
                     typename std::enable_if<not std::is_same<T, typename bagel::base_of<T>::type>::value
                                         and not std::is_const<T>::value
                                            >::type
                    > {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        typedef typename bagel::base_of<T>::type base;
        std::shared_ptr<base> u;
        load_impl<Archive, base>::load(ar, u);
        t = dynamic_pointer_cast<T>(u);
      }
    };

    // const classes
    template<class Archive, class T>
    struct load_impl<Archive, T,
                     typename std::enable_if<std::is_const<T>::value>::type
                    > {
      static void load(Archive& ar, std::shared_ptr<T>& t) {
        typedef typename std::remove_cv<T>::type non_const;
        std::shared_ptr<non_const> u;
        load_impl<Archive, non_const>::load(ar, u);
        t = u;
      }
    };


    template<class Archive, class T>
    inline void load(Archive& ar, std::shared_ptr<T>& t, const unsigned int) {
      load_impl<Archive, T>::load(ar, t);
    }


    // serialization of tuple
    template <size_t N>
    struct Serialize {
      template<class Archive, typename... Args>
      static void serialize(Archive& ar, std::tuple<Args...>& t, const unsigned int version) {
        ar & std::get<N-1>(t);
        Serialize<N-1>::serialize(ar, t, version);
      }
    };

    template<>
    struct Serialize<0LLU> {
      template<class Archive, typename... Args>
      static void serialize(Archive&, std::tuple<Args...>&, const unsigned int) { }
    };

    template<class Archive, typename... Args>
    void serialize(Archive& ar, std::tuple<Args...>& t, const unsigned int version) {
      Serialize<sizeof...(Args)>::serialize(ar, t, version);
    }

    // serizlization of std::array
    template<class Archive, typename T, size_t N>
    void serialize(Archive& ar, std::array<T,N>& t, const unsigned int) {
      for (size_t i = 0; i != N; ++i)
        ar & t[i];
    }

  }
}

#endif