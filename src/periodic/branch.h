//
// BAGEL - Parallel electron correlation program.
// Filename: branch.h
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_BRANCH_H
#define __SRC_PERIODIC_BRANCH_H

#include <src/molecule/shellpair.h>

namespace bagel {

class Branch {
  protected:
    int iws_;
    std::array<double, 3> centre_;
    std::vector<std::shared_ptr<const ShellPair>> sp_;

    double extent_;
    std::vector<std::shared_ptr<const ShellPair>> inter_;
    std::vector<std::shared_ptr<const Shellpair>> neigh_;

    bool is_neigh(std::shared_ptr<const Branch> b) const;
    void add_neigh(std::shared_ptr<const Branch> b);
    void add_inter(std::shared_ptr<const Branch> b);

  public:
    Branch(const int ws, const std::array<double, 3>& c, const std::vector<std::shared_ptr<const ShellPair>>& sp);
    ~Branch() { }

    const std::vector<std::shared_ptr<const ShellPair>>& neigh() const { return neigh_; }
    const std::vector<std::shared_ptr<const ShellPair>>& interaction_list() const { return inter_; }
    const std::vector<std::shared_ptr<const ShellPair>>& sp() const { return sp_; }
    std::shared_ptr<const ShellPair> sp(const int i) const { return sp_[i]; }
};

}
#endif