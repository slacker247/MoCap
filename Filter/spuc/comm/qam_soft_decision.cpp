#include <math.h>
#include <complex.h>
#include <qam_soft_decision.h>
namespace SPUC {
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
//!  
//!  \author Tony Kirke,  Copyright(c) 2001 
//! \ingroup comm
//! \brief  Soft-decision QAM Demapper for use with QAM_MOD class
//
//!  Returns an array of soft decision bits, array size is dependant on modulation type
//!  <p>This de-maps the QAM used in 802.11A
void qam_data_demap(long rate_index, complex<long> data_in,
					long soft_decision_level, long* viterbi_input)
{
  switch (rate_index) {
  case 1: (qpsk_soft_decision(data_in, viterbi_input));  break;
  case 2: (qam16_soft_decision(data_in, soft_decision_level, viterbi_input));  break;
  case 3: (qam64_soft_decision(data_in, soft_decision_level, viterbi_input));  break;
  default: (bpsk_soft_decision(data_in, viterbi_input));  break;
  }
}
//--------------------------------------------------------------------------
// BPSK data demapping
//--------------------------------------------------------------------------
void bpsk_soft_decision(complex<long> data_in, long* viterbi_input)
{
  viterbi_input[0] = (data_in.re); 
}
//--------------------------------------------------------------------------
// QPSK data demapping
//--------------------------------------------------------------------------
void qpsk_soft_decision(complex<long> data_in, long* viterbi_input)
{
  viterbi_input[0] = data_in.re;
  viterbi_input[1] = data_in.im;
}
//--------------------------------------------------------------------------
// 16 QAM data demapping
//--------------------------------------------------------------------------
void qam16_soft_decision(complex<long> data_in, long soft_decision_level, long* viterbi_input)
{
  viterbi_input[0] = data_in.re;
  if (data_in.re < 0) viterbi_input[1] = soft_decision_level + data_in.re;
  else 		          viterbi_input[1] = soft_decision_level - data_in.re;
  viterbi_input[2] = data_in.im;
  if (data_in.im < 0) viterbi_input[3] = soft_decision_level + data_in.im;
  else 		          viterbi_input[3] = soft_decision_level - data_in.im;
}
//--------------------------------------------------------------------------
// 64 QAM data demapping
//--------------------------------------------------------------------------
void qam64_soft_decision(complex<long> data_in, long soft_decision_level, long* viterbi_input)
{
  // real part
  viterbi_input[0] = data_in.re;
  if (data_in.re < 0) viterbi_input[1] = 2*soft_decision_level + data_in.re;
  else 		          viterbi_input[1] = 2*soft_decision_level - data_in.re;
  if (abs(data_in.re) < 2*soft_decision_level) {
	  if (data_in.re < 0) viterbi_input[2] = -soft_decision_level - data_in.re;
	  else 		          viterbi_input[2] = -soft_decision_level + data_in.re;
  } else {
	  if (data_in.re < 0) viterbi_input[2] = 3*soft_decision_level + data_in.re;
	  else 		          viterbi_input[2] = 3*soft_decision_level - data_in.re;
  }
  // imaginary part
  viterbi_input[3] = data_in.im;
  if (data_in.im < 0) viterbi_input[4] = 2*soft_decision_level + data_in.im;
  else 		          viterbi_input[4] = 2*soft_decision_level - data_in.im;
  if (abs(data_in.im) < 2*soft_decision_level) {
	  if (data_in.im < 0) viterbi_input[5] = -soft_decision_level - data_in.im;
	  else 		          viterbi_input[5] = -soft_decision_level + data_in.im;
  } else {
	  if (data_in.im < 0) viterbi_input[5] = 3*soft_decision_level + data_in.im;
	  else 		          viterbi_input[5] = 3*soft_decision_level - data_in.im;
  }
}
}
