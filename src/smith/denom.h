//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: denom.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SRC_SMITH_DENOM_H
#define __SRC_SMITH_DENOM_H

#include <src/wfn/rdm.h>
#include <src/util/math/zmatrix.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class Denom {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

  protected:
    std::shared_ptr<const MatType> fock_;
    const double thresh_;

    std::shared_ptr<MatType> shalf_x_;
    std::shared_ptr<MatType> shalf_h_;
    std::shared_ptr<MatType> shalf_xx_;
    std::shared_ptr<MatType> shalf_hh_;
    std::shared_ptr<MatType> shalf_xh_;
    std::shared_ptr<MatType> shalf_xhh_;
    std::shared_ptr<MatType> shalf_xxh_;

    std::shared_ptr<MatType> work_x_;
    std::shared_ptr<MatType> work_h_;
    std::shared_ptr<MatType> work_xx_;
    std::shared_ptr<MatType> work_hh_;
    std::shared_ptr<MatType> work_xh_;
    std::shared_ptr<MatType> work_xhh_;
    std::shared_ptr<MatType> work_xxh_;

    VectorB denom_x_;
    VectorB denom_h_;
    VectorB denom_xx_;
    VectorB denom_hh_;
    VectorB denom_xh_;
    VectorB denom_xhh_;
    VectorB denom_xxh_;

    // init functions
    void init_x_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                       std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_h_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                       std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xx_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_hh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xhh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                         std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xxh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                         std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);

  public:
    Denom(std::shared_ptr<const MatType> fock, const int nstates, const double thresh_overlap);

    // add RDMs (using fock-multiplied 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    // diagonalize and set to shalf and denom
    void compute();

    std::shared_ptr<const MatType> shalf_x() const { return shalf_x_; }
    std::shared_ptr<const MatType> shalf_h() const { return shalf_h_; }
    std::shared_ptr<const MatType> shalf_xx() const { return shalf_xx_; }
    std::shared_ptr<const MatType> shalf_hh() const { return shalf_hh_; }
    std::shared_ptr<const MatType> shalf_xh() const { return shalf_xh_; }
    std::shared_ptr<const MatType> shalf_xhh() const { return shalf_xhh_; }
    std::shared_ptr<const MatType> shalf_xxh() const { return shalf_xxh_; }

    const double& denom_x(const size_t i) const { return denom_x_(i); }
    const double& denom_h(const size_t i) const { return denom_h_(i); }
    const double& denom_xx(const size_t i) const { return denom_xx_(i); }
    const double& denom_hh(const size_t i) const { return denom_hh_(i); }
    const double& denom_xh(const size_t i) const { return denom_xh_(i); }
    const double& denom_xhh(const size_t i) const { return denom_xhh_(i); }
    const double& denom_xxh(const size_t i) const { return denom_xxh_(i); }

};

extern template class Denom<double>;
extern template class Denom<std::complex<double>>;

}
}

#endif
