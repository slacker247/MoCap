//! 
//! author="Tony Kirke" *
//! Copyright(c) 2001 Tony Kirke
// SPUC - Signal processing using C++ - A DSP library
/* 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
#include <spuc.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "qam_mod.h"
using namespace std;
using namespace SPUC;
//--------------------------------------------------------------------------
complex<long> qam_mod::data_map(long rate, long data)
{
  complex<long> data_out;
  switch (rate) {
  case 0: data_out = bpsk_map(data); break;
  case 1: data_out = qpsk_map(data); break;
  case 2: data_out = qam16_map(data); break;
  case 3: data_out = qam64_map(data); break;
  }
  return(data_out);
}
//--------------------------------------------------------------------------
// BPSK data mapping
//--------------------------------------------------------------------------
complex<long> qam_mod::bpsk_map(long data_in)
{
  if (data_in) return(complex<long>(1,0));
  else return(complex<long>(-1,0));
}
//--------------------------------------------------------------------------
// QPSK data mapping
//--------------------------------------------------------------------------
complex<long> qam_mod::qpsk_map(long phase_in)
{
  complex<long> data_out;
  switch (phase_in) {
  case 0: data_out = complex<long>(-1,-1); break;
  case 1: data_out = complex<long>(-1,1); break;
  case 2: data_out = complex<long>(1,-1); break;
  case 3: data_out = complex<long>(1,1); break;
  }
  return(data_out);
}
//--------------------------------------------------------------------------
// 16 QAM data mapping
//--------------------------------------------------------------------------
complex<long> qam_mod::qam16_map(long data_in)
{
  complex<long> data_out;
  switch (data_in) {
    case 0: data_out = complex<long>(-3,-3); break;
    case 1: data_out = complex<long>(-3,-1); break;
    case 2: data_out = complex<long>(-3,3); break;
    case 3: data_out = complex<long>(-3,1); break;
    case 4: data_out = complex<long>(-1,-3); break;
    case 5: data_out = complex<long>(-1,-1); break;
    case 6: data_out = complex<long>(-1,3); break;
    case 7: data_out = complex<long>(-1,1); break;
    case 8: data_out = complex<long>(3,-3); break;
    case 9: data_out = complex<long>(3,-1); break;
    case 10: data_out = complex<long>(3,3); break;
    case 11: data_out = complex<long>(3,1); break;
    case 12: data_out = complex<long>(1,-3); break;
    case 13: data_out = complex<long>(1,-1); break;
    case 14: data_out = complex<long>(1,3); break;
    case 15: data_out = complex<long>(1,1); break;
  }
  return(data_out);
}
//--------------------------------------------------------------------------
// 64 QAM data mapping
//--------------------------------------------------------------------------
complex<long> qam_mod::qam64_map(long data_in)
{
  long re,im;
  re = data_in >> 3;
  im = data_in & 0x7;

  switch (re) {
  case 0: re = -7; break;
  case 1: re = -5; break;
  case 2: re = -1; break;
  case 3: re = -3; break;
  case 4: re = 7; break;
  case 5: re = 5; break;
  case 6: re = 1; break;
  case 7: re = 3; break;
  default : 
	cout << "Illegal code re " << re << "\n";
	break;
  }
  switch (im) {
  case 0: im = -7; break;
  case 1: im = -5; break;
  case 2: im = -1; break;
  case 3: im = -3; break;
  case 4: im = 7; break;
  case 5: im = 5; break;
  case 6: im = 1; break;
  case 7: im = 3; break;
  default : 
	cout << "Illegal code im " << im << "\n";
	break;
  }
  return(complex<long>(re,im));
}
