//
// BAGEL - Parallel electron correlation program.
// Filename: reldffullt.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_DF_RELDFFULLT_H
#define __SRC_DF_RELDFFULLT_H

#include <src/df/dfdistt.h>
#include <src/df/reldffull.h>

namespace bagel {

class RelDFFullT {
  protected:
    std::array<std::shared_ptr<DFDistT>,2> dffull_;

  public:
    RelDFFullT(std::shared_ptr<const RelDFFull> full, std::shared_ptr<const StaticDist> dist = nullptr);

    const double* data(const int i) const { return dffull_.at(i)->data(); }

    // to make this class compatible with DFDistT
    int naux()   const { assert(dffull_[0]->naux() == dffull_[1]->naux());     return dffull_[0]->naux(); }
    int bstart() const { assert(dffull_[0]->bstart() == dffull_[1]->bstart()); return dffull_[0]->bstart(); }
    void discard_df();

    std::vector<std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>> get_slice(const int start, const int end) const;
    std::shared_ptr<ZMatrix> replicate() const;

    int locate(const size_t i, const size_t n) const { assert(dffull_[0]->locate(i,n) == dffull_[1]->locate(i,n)); return dffull_[0]->locate(i,n); }

    static int nblocks() { return 1; }
};

}

#endif