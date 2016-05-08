// 
//! Modified by Tony Kirke, Feb 1,2003
// SPUC - Signal processing using C++ - A DSP library
/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2001 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*! 
  \file 
  \brief Definition of a turbo encoder/decoder class
  \author Pål Frenger

  1.11 

  2002/11/05 00:40:02 
*/

#ifndef __turbo_h
#define __turbo_h

#include <binary.h>
#include <vector.h>
#include <matrix.h>
#include "rec_syst_conv_code.h"
#include "interleave.h"
namespace SPUC {
#ifndef DOXYGEN_SHOULD_SKIP_THIS
/*! 
  \brief Turbo encoder/decoder Class
  \ingroup fec

  To set up the turbo encoder used in e.g. WCDMA the following code can be used (assuming a code block size of 320 bits):
  Turbo_Codec turbo;
  ivec gen(2);
  gen(0) = 013; gen(1) = 015;
  int constraint_length = 4;
  ivec interleaver_sequence = wcdma_turbo_interleaver_sequence( 320 );
  turbo.set_parameters(gen, gen, constraint_length, interleaver_sequence);
*/
class Turbo_Codec {
 public:
  
  //! Class constructor
  Turbo_Codec(void) {}
  
  /*! 
    \brief Set parameters for the turbo encoder/decoder
    
    \param gen1 A vector with \a n1 elements containing the generator polynomials for the first constitiuent encoder 
    \param gen2 A vector with \a n2 elements containing the generator polynomials for the second constitiuent encoder 
    \param constraint_length The constraint length of the two constitiuent encoders 
    \param interleaver_sequence An ivec defining the internal turbo interleaver. 
    \param in_iterations The number decoding iterations. Default value is 8.
    \param in_metric Determines the decoder metric. May equal "MAP", LOGMAP", or "LOGMAX". Default value is "LOGMAX".
  */
  void set_parameters(ivec gen1, ivec gen2, int constraint_length, const ivec &interleaver_sequence, 
		      int in_iterations=8, string in_metric="LOGMAX");

  /*!
    \brief Set a new internal interleaver sequence for the turbo encoder/decoder

    By changing the interleaver sequence it is possible to change the code word size of the turbo codec. Note that you
    still must use the \a set_parameters function first to set the other parameters of the turbo codec.
   */
  void set_interleaver(const ivec &interleaver_sequence);
  
  /*! 
    \brief Set parameters for decoding on an AWGN channel

    \param in_Ec The received energy per channel symbol (i.e. per soft input bit)
    \param in_N0 The single sided spectral density of the AWGN noise
  */
  void set_awgn_channel_parameters(double in_Ec, double in_N0);
  
  /*! 
    \brief Set scaling factor for decoding on e.g. Rayleigh fading channels

    Setting the correct value of the channel reliability function is important for MAP decoder algorithms. However, if
    the Log-MAX decoding algorithm is used, then the value of \a Lc is not important. Assuming that the received soft 
    values \f$r_k\f$ from the channel equal

    \f[ r_k = h_k c_k + w_k \f]

    where \f$h_k\f$ is the (complex valued) channel gain, \f$c_k\f$ is the transmitted symbol equal to 
    \f$\{-\sqrt{E_c},+\sqrt{E_c}\}\f$, and \f$w_k\f$ is the AWGN (complex valued) noise with total variance
    of the real plus imaginary part equal to \f$N_0\f$. The input to the turbo decoder shall then be

    \f[ z_k = \hat{h}_k^{*} r_k \f]

    where \f$\hat{h}_k^{*}\f$ is the conjugate of the channel estimate. Assuming that the channel estimate is perfect,
    the channel reliability factor shall be set to

    \f[ L_c = 4\sqrt{E_c} / {N_0} \f]

    \param in_Lc the channel reliability factor of the channel.
  */
  void set_scaling_factor(double in_Lc);

  /*!
    \brief Encoder function

    This function can encode several consequtive coding blocks. The output is organized as follows:

    \f[ s(1), p_{1,1}(1), p_{1,2}(1), \ldots , p_{1,n_1}(1), p_{2,1}(1), p_{2,2}(1), \ldots , p_{2,n_2}(1), s(2), \ldots \f]

    In the above expression \f$s(n)\f$ is the n-th systematic bit and \f$p_{l,k}(n)\f$ is the n-th bit from the k-th encoder polynom
    of the l-th constituent encoder. A tail of both systematic and parity bits is added after each coding block to force both encoder
    to the zero state. The tail of each encoder is structured as follows (using encoder one as an example):

    \f[ t_1(1), pt_{1,1}(1), pt_{1,2}(1), \ldots , pt_{1,n_1}(1), \ldots pt_{1,n_1}(m) \f]

    The tailbits from the first encoder are placed before the tailbits from the second encoder.

    \param input The input bits the the encoder. Must contain an integer number of code blocks
    \param output The encoded bits including two tail, one from each constituent encoder, after each coding block.
  */
  void encode(const bvec &input, bvec &output);

  /*! 
    \brief Decoder function

    This function can decode several consequtive coding blocks that were encoded with the encode member function
  */
  void decode(const vec &received_signal, bvec &decoded_bits, const bvec &true_bits="0");

  /*! 
    \brief Encode a single block

    This function is usefull if rate matching is to be aplied after the turbo encoder. The size of \a in1 and \a in2
    equals the size of \a input plus the tail bits. Tailbits are appended ate the end of \a in1 and \a in2. The number
    of rows in \a parity1 and \a parity2 equals the lengths of \a in1 and \a in2 respectively. The number of columns of
    \a parity1 and \a parity2 equals the number of parity bits from the first and the second constituent encoders
    respectively.

    \param input The input bits the the encoder. Must contain a single code block
    \param in1 The input bits to the first constituent encoder with a tail added at the end
    \param in2 The input bits to the second constituent encoder (i.e. the interleaved data bits) with a tail added at the end
    \param parity1 The parity bits from the first constituent encoder (including parity bits for the first tail)
    \param parity2 The parity bits from the second constituent encoder (including parity bits for the second tail)
  */
  void encode_block(const bvec &input, bvec &in1, bvec &in2, bmat &parity1, bmat &parity2);

  /*! 
    \brief Decode a single block

    This function can decode a single coding blocks that was encoded with the encode_block member function. In order to speed
    up the decoding process it is possible to input the correct data bits. If this information is provided the decoder can
    stop iterating as soon as the decoded bits match the correct data bits. This simulation trick can greatly speed up the
    simulation time for high SNR cases when only a few iterations are required. If errors still exist after the maximum
    number of iterations have been performed, the decoding process stops.

    The matrix \a decoded_bits_i contains the result from decoding iteration \a i on row \a i. Even though both \a rec_syst1 
    and \a rec_syst2 are given as inputs, the systematic bits in \a rec_syst2 will in most cases be punctured and only the 
    tailbits at the end of the vector \a rec_syst2 will have values different from zero. 

    \note This decoding function assumes that the input is scaled as +-2*SNR + noise. This means that the channel 
    reliability factor \a Lc must be equal to 1.0. No additional scaling is performed by this function.

    \param rec_syst1 The received input bits to the first constituent decoder with a tail added at the end
    \param rec_syst2 The received input bits to the second constituent decoder with a tail added at the end
    \param rec_parity1 The received parity bits for the first constituent decoder (including parity bits for the first tail)
    \param rec_parity2 The received parity bits for the second constituent decoder (including parity bits for the second tail)
    \param decoded_bits_i Contains the result from decoding iteration \a i on row \a i. 
    \param true_bits Optional input parameter. If given, the iterations will stop as soon as the decoded bits match the true bits.
  */
  void decode_block(const vec &rec_syst1, const vec &rec_syst2, const mat &rec_parity1, const mat &rec_parity2, 
		    bmat &decoded_bits_i, const bvec &true_bits="0");

private:   

  /*! 
    \brief Special decoder function for \a n = 2
  */
  void decode_n2(const vec &received_signal, bvec &decoded_bits);

  //Scalars:
  long interleaver_size;
  long Ncoded, Nuncoded;
  int m_tail, n1, n2, n_tot, iterations;
  double Ec, N0, Lc, R;
  string metric;

  //Classes:
  Rec_Syst_Conv_Code rscc1, rscc2;
  Sequence_Interleaver<bin> bit_interleaver;
  Sequence_Interleaver<double> float_interleaver;

};


/*!
  \relates Turbo_Codec
  \brief Generates the interleaver sequence for the internal turbo encoder interleaver used in WCDMA
*/
ivec wcdma_turbo_interleaver_sequence(int interleaver_size);
#endif
}
#endif // __turbo_h
