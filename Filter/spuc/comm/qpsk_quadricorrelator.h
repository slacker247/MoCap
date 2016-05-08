namespace SPUC {
/*
 * SPUC - Signal processing using C++ - A DSP library
 * 
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
//!	  \ingroup comm
//! \brief  QPSK Quadricorrelator frequency discriminator
//!  \author Tony Kirke,  Copyright(c) 2001 
//!  
//
//!   Hard_data - Hard decision (complex) data for current symbol
//!   Prev - Previous symbol (prior to hard decision)
//! \image html qpsk_quadricorrelator.gif
//! \image latex qpsk_quadricorrelator.eps
template <class T> T qpsk_quadricorrelator(complex<T> hard_data, complex<T> prev) 
{
	return(hard_data.re*prev.im-hard_data.im*prev.re);
}
} // namespace SPUC 
