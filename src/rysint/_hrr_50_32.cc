//
// BAGEL - Parallel electron correlation program.
// Filename: _hrr_50_32.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/rysint/hrrlist.h>
#include <array>
#include <algorithm>

using namespace std;
using namespace bagel;

void HRRList::perform_HRR_50_32(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 46];
    double* current_out = &data_out[c * 60];
   {
     //current index a: xxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[10];
      const double ay_0 = current_data[11];
      const double az_0 = current_data[15];
      const double axx_0 = current_data[25];
      const double axy_0 = current_data[26];
      const double ayy_0 = current_data[27];
      const double axz_0 = current_data[31];
      const double ayz_0 = current_data[32];
      const double azz_0 = current_data[36];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[0] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[1] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[2] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[3] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[4] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[5] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[11];
      const double ay_0 = current_data[12];
      const double az_0 = current_data[16];
      const double axx_0 = current_data[26];
      const double axy_0 = current_data[27];
      const double ayy_0 = current_data[28];
      const double axz_0 = current_data[32];
      const double ayz_0 = current_data[33];
      const double azz_0 = current_data[37];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[6] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[7] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[8] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[9] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[10] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[11] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[12];
      const double ay_0 = current_data[13];
      const double az_0 = current_data[17];
      const double axx_0 = current_data[27];
      const double axy_0 = current_data[28];
      const double ayy_0 = current_data[29];
      const double axz_0 = current_data[33];
      const double ayz_0 = current_data[34];
      const double azz_0 = current_data[38];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[12] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[13] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[14] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[15] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[16] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[17] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[13];
      const double ay_0 = current_data[14];
      const double az_0 = current_data[18];
      const double axx_0 = current_data[28];
      const double axy_0 = current_data[29];
      const double ayy_0 = current_data[30];
      const double axz_0 = current_data[34];
      const double ayz_0 = current_data[35];
      const double azz_0 = current_data[39];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[18] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[19] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[20] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[21] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[22] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[23] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxz
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[15];
      const double ay_0 = current_data[16];
      const double az_0 = current_data[19];
      const double axx_0 = current_data[31];
      const double axy_0 = current_data[32];
      const double ayy_0 = current_data[33];
      const double axz_0 = current_data[36];
      const double ayz_0 = current_data[37];
      const double azz_0 = current_data[40];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[24] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[25] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[26] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[27] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[28] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[29] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyz
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[16];
      const double ay_0 = current_data[17];
      const double az_0 = current_data[20];
      const double axx_0 = current_data[32];
      const double axy_0 = current_data[33];
      const double ayy_0 = current_data[34];
      const double axz_0 = current_data[37];
      const double ayz_0 = current_data[38];
      const double azz_0 = current_data[41];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[30] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[31] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[32] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[33] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[34] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[35] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyz
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[17];
      const double ay_0 = current_data[18];
      const double az_0 = current_data[21];
      const double axx_0 = current_data[33];
      const double axy_0 = current_data[34];
      const double ayy_0 = current_data[35];
      const double axz_0 = current_data[38];
      const double ayz_0 = current_data[39];
      const double azz_0 = current_data[42];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[36] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[37] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[38] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[39] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[40] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[41] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xzz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[19];
      const double ay_0 = current_data[20];
      const double az_0 = current_data[22];
      const double axx_0 = current_data[36];
      const double axy_0 = current_data[37];
      const double ayy_0 = current_data[38];
      const double axz_0 = current_data[40];
      const double ayz_0 = current_data[41];
      const double azz_0 = current_data[43];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[42] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[43] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[44] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[45] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[46] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[47] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yzz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[20];
      const double ay_0 = current_data[21];
      const double az_0 = current_data[23];
      const double axx_0 = current_data[37];
      const double axy_0 = current_data[38];
      const double ayy_0 = current_data[39];
      const double axz_0 = current_data[41];
      const double ayz_0 = current_data[42];
      const double azz_0 = current_data[44];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[48] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[49] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[50] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[51] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[52] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[53] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: zzz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[22];
      const double ay_0 = current_data[23];
      const double az_0 = current_data[24];
      const double axx_0 = current_data[40];
      const double axy_0 = current_data[41];
      const double ayy_0 = current_data[42];
      const double axz_0 = current_data[43];
      const double ayz_0 = current_data[44];
      const double azz_0 = current_data[45];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[54] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[55] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[56] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[57] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[58] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[59] = az_z + AB[2] * a0_z; // a0_zz

    }
  }
}

