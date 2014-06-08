//
// BAGEL - Parallel electron correlation program.
// Filename: dfdistt_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_DF_DFDISTT_LONDON_H
#define __SRC_DF_DFDISTT_LONDON_H

#include <src/df/df_london.h>

namespace bagel {

/*
    DFDistT_London is a 3-index DF integral object, which is distributed by the second and third indices.
*/

class DFDistT_London {
  protected:
    std::vector<std::shared_ptr<ZMatrix>> data_;

    // first dimension is naux_ (global)
    const size_t naux_;
    // original dimensions
    const size_t nindex1_;
    const size_t nindex2_;

    // distribution information
    std::shared_ptr<const StaticDist> dist_;

    // second and third dimension
    size_t bstart_;
    size_t bsize_;

    std::shared_ptr<const ParallelDF_London> df_;

  public:
    // CAUTION this constructor should be called **COLLECTIVELY**!! Otherwise the program hangs.
    DFDistT_London(std::shared_ptr<const ParallelDF_London> in, std::shared_ptr<const StaticDist> dist = nullptr);

    DFDistT_London(const size_t naux, std::shared_ptr<const StaticDist> dist, const size_t n1, const size_t n2,
            const std::shared_ptr<const ParallelDF_London>);

    std::shared_ptr<DFDistT_London> clone() const;
    std::shared_ptr<DFDistT_London> apply_J(std::shared_ptr<const ZMatrix> d) const;
    std::shared_ptr<DFDistT_London> apply_J() const { return apply_J(df_->data2()); }
    std::vector<std::shared_ptr<ZMatrix>> form_aux_2index(std::shared_ptr<const DFDistT_London> o, const double a) const;

    void get_paralleldf(std::shared_ptr<ParallelDF_London>) const;

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }

    int bsize() const { return bsize_; }
    int bstart() const { return bstart_; }
    int nblocks() const { return data_.size(); }
    const std::complex<double>* data() const { assert(data_.size() == 1); return data(0); }
    const std::complex<double>* data(const int i) const { return data_[i]->data(); }

    // returns the process that has the data
    int locate(const size_t, const size_t n) const { return std::get<0>(dist_->locate(n)); }

    std::vector<std::shared_ptr<ZMatrix>> get_slice(const int start, const int end) const;

    std::shared_ptr<const ParallelDF_London> df() const { return df_; }
    void discard_df() { df_.reset(); }

    std::shared_ptr<ZMatrix> replicate(const int i = 0) const;
};

}

#endif