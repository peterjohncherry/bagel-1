//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_tasks2.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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


#include <src/smith/CASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

void Task50::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I57
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  target_ += ddot_(c1.size()*a2.size()*c3.size()*a4.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task51::Task_local::compute() {
  target_ = 0.0;
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I57
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a2, c3, a4);
}

void Task52::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: Gamma0
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I59
  std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x1)]);
  sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x1.size());
  target_ += ddot_(x0.size()*x1.size(), i0data_sorted, 1, i1data_sorted, 1);
}

void Task53::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I59
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I60
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task54::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I60
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task55::Task_local::compute() {
  target_ = 0.0;
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: I59
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x1), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I63
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x1.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x0, x1);
}

void Task56::Task_local::compute() {
  target_ = 0.0;
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I63
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task58::Task_local::compute() {
  const Index x0 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, x1);
  {
    // tensor label: I64
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x1.size(), x0.size());
  }
  out()->put_block(odata, x0, x1);
}

void Task59::Task_local::compute() {
  const Index x1 = b(0);
  const Index x0 = b(1);
  // tensor label: I64
  std::unique_ptr<double[]> odata = out()->move_block(x1, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x1, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x1, x0), 0.0);
  // tensor label: Gamma0
  std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
  sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
  // tensor label: I65
  std::unique_ptr<double[]> i1data = in(1)->get_block();
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size()]);
  sort_indices<0,1,1,1>(i1data, i1data_sorted);
  dgemm_("T", "N", x1.size()*x0.size(), 1, 1,
         1.0, i0data_sorted, 1, i1data_sorted, 1,
         1.0, odata_sorted, x1.size()*x0.size());
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x1.size(), x0.size());
  out()->put_block(odata, x1, x0);
}

void Task60::Task_local::compute() {
  // tensor label: I65
  std::unique_ptr<double[]> odata = out()->move_block();
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index a2 = b(0);
  const Index c3 = b(1);
  const Index a4 = b(2);
  const Index c1 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c3, a2);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a4, c3, a2)]);
  sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a4.size(), c3.size(), a2.size());
  // tensor label: I66
  std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
  sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
  odata_sorted[0] += ddot_(a4.size()*c3.size()*a2.size()*c1.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->put_block(odata);
}

void Task61::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I66
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task62::Task_local::compute() {
  // tensor label: I65
  std::unique_ptr<double[]> odata = out()->move_block();
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size()]);
  std::fill_n(odata_sorted.get(), out()->get_size(), 0.0);
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: t2
  std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  // tensor label: I69
  std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a2, c3, a4);
  std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a2, c3, a4)]);
  sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
  odata_sorted[0] += ddot_(c1.size()*a2.size()*c3.size()*a4.size(), i0data_sorted, 1, i1data_sorted, 1);
  sort_indices<1,1,1,1>(odata_sorted, odata);
  out()->put_block(odata);
}

void Task63::Task_local::compute() {
  const Index c1 = b(0);
  const Index a2 = b(1);
  const Index c3 = b(2);
  const Index a4 = b(3);
  // tensor label: I69
  std::unique_ptr<double[]> odata = out()->move_block(c1, a2, c3, a4);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, c1, a2, c3, a4);
}

void Task64::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  {
    // tensor label: I70
    std::unique_ptr<double[]> i0data = in(0)->get_block(c5, c3);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c5.size(), c3.size());
  }
  out()->put_block(odata, c5, c3);
}

void Task65::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: I70
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c1 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a4)]);
        sort_indices<3,1,0,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a4.size());
        // tensor label: I71
        std::unique_ptr<double[]> i1data = in(1)->get_block(c1, a4, c5, a2);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c1, a4, c5, a2)]);
        sort_indices<1,3,0,2,0,1,1,1>(i1data, i1data_sorted, c1.size(), a4.size(), c5.size(), a2.size());
        dgemm_("T", "N", c3.size(), c5.size(), c1.size()*a4.size()*a2.size(),
               1.0, i0data_sorted, c1.size()*a4.size()*a2.size(), i1data_sorted, c1.size()*a4.size()*a2.size(),
               1.0, odata_sorted, c3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c3.size(), c5.size());
  out()->put_block(odata, c5, c3);
}

void Task66::Task_local::compute() {
  const Index c1 = b(0);
  const Index a4 = b(1);
  const Index c5 = b(2);
  const Index a2 = b(3);
  // tensor label: I71
  std::unique_ptr<double[]> odata = out()->move_block(c1, a4, c5, a2);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a4, c5, a2);
    sort_indices<0,1,2,3,1,1,2,1>(i0data, odata, c1.size(), a4.size(), c5.size(), a2.size());
  }
  out()->put_block(odata, c1, a4, c5, a2);
}

void Task67::Task_local::compute() {
  const Index c5 = b(0);
  const Index c3 = b(1);
  // tensor label: I70
  std::unique_ptr<double[]> odata = out()->move_block(c5, c3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c5, c3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c5, c3), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c5, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c5, a4)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c5.size(), a4.size());
        // tensor label: I73
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", c5.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, c5.size());
      }
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, c5.size(), c3.size());
  out()->put_block(odata, c5, c3);
}

void Task68::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I73
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-4,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task69::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
  {
    // tensor label: I74
    std::unique_ptr<double[]> i0data = in(0)->get_block(a4, a5);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a4.size(), a5.size());
  }
  out()->put_block(odata, a4, a5);
}

void Task70::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& c3 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a5, c3, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a5, c3, a2)]);
        sort_indices<0,2,3,1,0,1,1,1>(i0data, i0data_sorted, c1.size(), a5.size(), c3.size(), a2.size());
        // tensor label: I75
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<3,1,2,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a5.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a4, a5);
}

void Task71::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I75
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-2,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task72::Task_local::compute() {
  const Index a4 = b(0);
  const Index a5 = b(1);
  // tensor label: I74
  std::unique_ptr<double[]> odata = out()->move_block(a4, a5);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a4, a5)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a4, a5), 0.0);
  for (auto& c1 : *range_[0]) {
    for (auto& a2 : *range_[2]) {
      for (auto& c3 : *range_[0]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a5);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c1, a2, c3, a5)]);
        sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, c1.size(), a2.size(), c3.size(), a5.size());
        // tensor label: I77
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<3,2,1,0,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", a5.size(), a4.size(), c3.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, c3.size()*a2.size()*c1.size(), i1data_sorted, c3.size()*a2.size()*c1.size(),
               1.0, odata_sorted, a5.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, a5.size(), a4.size());
  out()->put_block(odata, a4, a5);
}

void Task73::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I77
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,4,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task74::Task_local::compute() {
  const Index x0 = b(0);
  const Index c3 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x0, c3);
  {
    // tensor label: I78
    std::unique_ptr<double[]> i0data = in(0)->get_block(c3, x0);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c3.size(), x0.size());
  }
  out()->put_block(odata, x0, c3);
}

void Task75::Task_local::compute() {
  const Index c3 = b(0);
  const Index x0 = b(1);
  // tensor label: I78
  std::unique_ptr<double[]> odata = out()->move_block(c3, x0);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x0)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x0), 0.0);
  for (auto& x1 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<0,1,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I79
    std::unique_ptr<double[]> i1data = in(1)->get_block(c3, x1);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(c3, x1)]);
    sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, c3.size(), x1.size());
    dgemm_("T", "N", x0.size(), c3.size(), x1.size(),
           1.0, i0data_sorted, x1.size(), i1data_sorted, x1.size(),
           1.0, odata_sorted, x0.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x0.size(), c3.size());
  out()->put_block(odata, c3, x0);
}

void Task76::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& a4 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a2 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a4, c1, a2);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a4, c1, a2)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a4.size(), c1.size(), a2.size());
        // tensor label: I80
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<0,3,2,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x1.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task77::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I80
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task78::Task_local::compute() {
  const Index c3 = b(0);
  const Index x1 = b(1);
  // tensor label: I79
  std::unique_ptr<double[]> odata = out()->move_block(c3, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c3, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c3, x1), 0.0);
  for (auto& a2 : *range_[2]) {
    for (auto& c1 : *range_[0]) {
      for (auto& a4 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a2, c1, a4);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a2, c1, a4)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), a2.size(), c1.size(), a4.size());
        // tensor label: I83
        std::unique_ptr<double[]> i1data = in(1)->get_block(a4, c3, a2, c1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a4, c3, a2, c1)]);
        sort_indices<2,3,0,1,0,1,1,1>(i1data, i1data_sorted, a4.size(), c3.size(), a2.size(), c1.size());
        dgemm_("T", "N", x1.size(), c3.size(), a4.size()*a2.size()*c1.size(),
               1.0, i0data_sorted, a4.size()*a2.size()*c1.size(), i1data_sorted, a4.size()*a2.size()*c1.size(),
               1.0, odata_sorted, x1.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c3.size());
  out()->put_block(odata, c3, x1);
}

void Task79::Task_local::compute() {
  const Index a4 = b(0);
  const Index c3 = b(1);
  const Index a2 = b(2);
  const Index c1 = b(3);
  // tensor label: I83
  std::unique_ptr<double[]> odata = out()->move_block(a4, c3, a2, c1);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(c1, a2, c3, a4);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, c1.size(), a2.size(), c3.size(), a4.size());
  }
  out()->put_block(odata, a4, c3, a2, c1);
}

void Task80::Task_local::compute() {
  const Index c4 = b(0);
  const Index x1 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c4, x1);
  {
    // tensor label: I84
    std::unique_ptr<double[]> i0data = in(0)->get_block(c4, x1);
    sort_indices<0,1,1,1,1,1>(i0data, odata, c4.size(), x1.size());
  }
  out()->put_block(odata, c4, x1);
}

void Task81::Task_local::compute() {
  const Index c4 = b(0);
  const Index x1 = b(1);
  // tensor label: I84
  std::unique_ptr<double[]> odata = out()->move_block(c4, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c4, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c4, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I85
    std::unique_ptr<double[]> i1data = in(1)->get_block(x0, c4);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, c4)]);
    sort_indices<0,1,0,1,1,1>(i1data, i1data_sorted, x0.size(), c4.size());
    dgemm_("T", "N", x1.size(), c4.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x1.size(), c4.size());
  out()->put_block(odata, c4, x1);
}

void Task82::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  // tensor label: I85
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a3.size(), c4.size(), a1.size());
        // tensor label: I86
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<1,0,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", c4.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size());
  out()->put_block(odata, x0, c4);
}

void Task83::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I86
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,1>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task84::Task_local::compute() {
  const Index x0 = b(0);
  const Index c4 = b(1);
  // tensor label: I85
  std::unique_ptr<double[]> odata = out()->move_block(x0, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, c4), 0.0);
  for (auto& c2 : *range_[0]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(c2, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(c2, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, c2.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I89
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<1,2,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", c4.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), x0.size());
  out()->put_block(odata, x0, c4);
}

void Task85::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I89
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task86::Task_local::compute() {
  const Index x1 = b(0);
  const Index x2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(x1, x2);
  {
    // tensor label: I90
    std::unique_ptr<double[]> i0data = in(0)->get_block(x2, x1);
    sort_indices<1,0,1,1,1,1>(i0data, odata, x2.size(), x1.size());
  }
  out()->put_block(odata, x1, x2);
}

void Task87::Task_local::compute() {
  const Index x2 = b(0);
  const Index x1 = b(1);
  // tensor label: I90
  std::unique_ptr<double[]> odata = out()->move_block(x2, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x2, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x2, x1), 0.0);
  for (auto& x3 : *range_[1]) {
    for (auto& x0 : *range_[1]) {
      // tensor label: Gamma26
      std::unique_ptr<double[]> i0data = in(0)->get_block(x3, x0, x2, x1);
      std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, x0, x2, x1)]);
      sort_indices<0,1,2,3,0,1,1,1>(i0data, i0data_sorted, x3.size(), x0.size(), x2.size(), x1.size());
      // tensor label: I91
      std::unique_ptr<double[]> i1data = in(1)->get_block(x0, x3);
      std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(x0, x3)]);
      sort_indices<1,0,0,1,1,1>(i1data, i1data_sorted, x0.size(), x3.size());
      dgemm_("T", "N", x2.size()*x1.size(), 1, x0.size()*x3.size(),
             1.0, i0data_sorted, x0.size()*x3.size(), i1data_sorted, x0.size()*x3.size(),
             1.0, odata_sorted, x2.size()*x1.size());
    }
  }
  sort_indices<0,1,1,1,1,1>(odata_sorted, odata, x2.size(), x1.size());
  out()->put_block(odata, x2, x1);
}

void Task88::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I91
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a3 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a3, c2, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a3, c2, a1)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a3.size(), c2.size(), a1.size());
        // tensor label: I92
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<0,1,2,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task89::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I92
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task90::Task_local::compute() {
  const Index x0 = b(0);
  const Index x3 = b(1);
  // tensor label: I91
  std::unique_ptr<double[]> odata = out()->move_block(x0, x3);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(x0, x3)]);
  std::fill_n(odata_sorted.get(), out()->get_size(x0, x3), 0.0);
  for (auto& a1 : *range_[2]) {
    for (auto& c2 : *range_[0]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x3, a1, c2, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x3, a1, c2, a3)]);
        sort_indices<1,2,3,0,0,1,1,1>(i0data, i0data_sorted, x3.size(), a1.size(), c2.size(), a3.size());
        // tensor label: I95
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
        sort_indices<2,1,0,3,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
        dgemm_("T", "N", x3.size(), x0.size(), a3.size()*c2.size()*a1.size(),
               1.0, i0data_sorted, a3.size()*c2.size()*a1.size(), i1data_sorted, a3.size()*c2.size()*a1.size(),
               1.0, odata_sorted, x3.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, x3.size(), x0.size());
  out()->put_block(odata, x0, x3);
}

void Task91::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I95
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task92::Task_local::compute() {
  const Index c4 = b(0);
  const Index c2 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(c4, c2);
  {
    // tensor label: I96
    std::unique_ptr<double[]> i0data = in(0)->get_block(c2, c4);
    sort_indices<1,0,1,1,1,1>(i0data, odata, c2.size(), c4.size());
  }
  out()->put_block(odata, c4, c2);
}

void Task93::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I96
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a3 : *range_[2]) {
      for (auto& a1 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a3, c4, a1);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a3, c4, a1)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a3.size(), c4.size(), a1.size());
        // tensor label: I97
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,0,2,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c2.size(), a3.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*a1.size()*x1.size(), i1data_sorted, a3.size()*a1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c2.size());
  out()->put_block(odata, c2, c4);
}

void Task94::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I97
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I98
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task95::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I98
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,1,4>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task96::Task_local::compute() {
  const Index c2 = b(0);
  const Index c4 = b(1);
  // tensor label: I96
  std::unique_ptr<double[]> odata = out()->move_block(c2, c4);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(c2, c4)]);
  std::fill_n(odata_sorted.get(), out()->get_size(c2, c4), 0.0);
  for (auto& x1 : *range_[1]) {
    for (auto& a1 : *range_[2]) {
      for (auto& a3 : *range_[2]) {
        // tensor label: t2
        std::unique_ptr<double[]> i0data = in(0)->get_block(x1, a1, c4, a3);
        std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, a1, c4, a3)]);
        sort_indices<0,1,3,2,0,1,1,1>(i0data, i0data_sorted, x1.size(), a1.size(), c4.size(), a3.size());
        // tensor label: I100
        std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x1);
        std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x1)]);
        sort_indices<3,2,0,1,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x1.size());
        dgemm_("T", "N", c4.size(), c2.size(), a3.size()*a1.size()*x1.size(),
               1.0, i0data_sorted, a3.size()*a1.size()*x1.size(), i1data_sorted, a3.size()*a1.size()*x1.size(),
               1.0, odata_sorted, c4.size());
      }
    }
  }
  sort_indices<1,0,1,1,1,1>(odata_sorted, odata, c4.size(), c2.size());
  out()->put_block(odata, c2, c4);
}

void Task97::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x1 = b(3);
  // tensor label: I100
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x1);
  std::unique_ptr<double[]> odata_sorted(new double[out()->get_size(a3, c2, a1, x1)]);
  std::fill_n(odata_sorted.get(), out()->get_size(a3, c2, a1, x1), 0.0);
  for (auto& x0 : *range_[1]) {
    // tensor label: Gamma0
    std::unique_ptr<double[]> i0data = in(0)->get_block(x1, x0);
    std::unique_ptr<double[]> i0data_sorted(new double[in(0)->get_size(x1, x0)]);
    sort_indices<1,0,0,1,1,1>(i0data, i0data_sorted, x1.size(), x0.size());
    // tensor label: I101
    std::unique_ptr<double[]> i1data = in(1)->get_block(a3, c2, a1, x0);
    std::unique_ptr<double[]> i1data_sorted(new double[in(1)->get_size(a3, c2, a1, x0)]);
    sort_indices<3,0,1,2,0,1,1,1>(i1data, i1data_sorted, a3.size(), c2.size(), a1.size(), x0.size());
    dgemm_("T", "N", x1.size(), a3.size()*c2.size()*a1.size(), x0.size(),
           1.0, i0data_sorted, x0.size(), i1data_sorted, x0.size(),
           1.0, odata_sorted, x1.size());
  }
  sort_indices<1,2,3,0,1,1,1,1>(odata_sorted, odata, x1.size(), a3.size(), c2.size(), a1.size());
  out()->put_block(odata, a3, c2, a1, x1);
}

void Task98::Task_local::compute() {
  const Index a3 = b(0);
  const Index c2 = b(1);
  const Index a1 = b(2);
  const Index x0 = b(3);
  // tensor label: I101
  std::unique_ptr<double[]> odata = out()->move_block(a3, c2, a1, x0);
  {
    // tensor label: t2
    std::unique_ptr<double[]> i0data = in(0)->get_block(x0, a1, c2, a3);
    sort_indices<3,2,1,0,1,1,-1,2>(i0data, odata, x0.size(), a1.size(), c2.size(), a3.size());
  }
  out()->put_block(odata, a3, c2, a1, x0);
}

void Task99::Task_local::compute() {
  const Index a1 = b(0);
  const Index a4 = b(1);
  // tensor label: r
  std::unique_ptr<double[]> odata = out()->move_block(a1, a4);
  {
    // tensor label: I102
    std::unique_ptr<double[]> i0data = in(0)->get_block(a1, a4);
    sort_indices<0,1,1,1,1,1>(i0data, odata, a1.size(), a4.size());
  }
  out()->put_block(odata, a1, a4);
}

