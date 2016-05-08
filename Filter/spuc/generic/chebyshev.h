// 
// author="Tony Kirke"
// Copyright(c) 1993-1998 Tony Kirke
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
#include <spuc.h>
#include <complex.h>
#include <iir_1st.h>
#include <iir_2nd.h>
namespace SPUC {
/*! 
  \addtogroup iir IIR filters
*/

/*!  \brief  Template Class for Chebyshev low pass iir filter
  \ingroup iir
*/
//! Chebyshev low pass filter.
//! Specify cut-off, order and passband ripple on construction
//! Filter is bilinear transformation of analog chebyshev
template <class Numeric> class chebyshev 
{
	private:
	long order;
	long odd;
	long n2;
	complex<double>* roots;
	iir_2nd< Numeric>* iir;
	iir_1st< Numeric>* iir_1;
	double gain;
	double* coeff;
	double epi;

	public:
	//! Constructor, fcd = cut-off (1=sampling rate)
	//! ord = Filter order
	//! ripple = passband ripple in dB
	chebyshev(double fcd, long ord=1, double ripple=3.0) {
		gain = 1;
		order = ord;
		epi = pow(10.0,(ripple/10.)) - 1.0;
		epi = pow(epi,(1./(1.0*ord)));
		double wca = tan(0.5*PI*fcd);
		//! wca - pre-warped angular frequency
	    n2 = (order+1)/2;
		odd = (order%2);
		roots = new complex<double>[n2];
		if (odd) iir_1 = new iir_1st< Numeric >;
		get_roots(wca,order,n2);
		bilinear(n2);
		iir = new iir_2nd< Numeric >[order/2];
		get_coeff(n2);
	}
	//! Destructor
	~chebyshev() {
	  delete [] roots;
	  if (odd) delete iir_1;
	  delete [] iir;
	}
	//! Reset history
	void reset() {
	  for (int j=odd;j<n2;j++) { 
		iir[j-odd].reset();
		if (odd) iir_1->reset();
	  }
	} 
	//! print coefficients
	void print() {
		int j;
		for (j=odd;j<n2;j++) { 
			iir[j-odd].print();
		}
		if (odd) iir_1->print();
	} 
	//! Clock in sample and get output.
	Numeric clock(Numeric in) {
		Numeric tmp = in;
		for (int i=odd;i<n2;i++) {
			tmp = iir[i-odd].clock(tmp);
		}
		if (odd) tmp = iir_1->clock(tmp);
		return(gain*tmp);
	}
	private:
	//! Calculate roots (chebyshev)
	void get_roots(double wp, long n, long n2) {
       long l = 0;
       if (n%2 == 0) l = 1;                                                 
	   double arg;
	   double x = 1/epi;
	   double asinh = log(x + sqrt(1.0+x*x));
	   double v0 = asinh/(double(n));
	   double sm = sinh(v0);
	   double cm = cosh(v0);
	   for (int j=0;j<n2;j++) {
		 arg = -0.5*PI*l/float(n);
		 roots[j] = wp*complex<double>(-sm*cos(arg),cm*sin(arg));
		 l += 2;
	  }
	}
	//! Do bilinear transformation
	void bilinear(long n2) {
      double td;
	  int j;
	  const double a = 1.;
      if (odd) {
		// For s-domain zero (or pole) at infinity we
		// get z-domain zero (or pole) at z = -1
		double tmp = (roots[0].real()-a)/(a+roots[0].real());
		iir_1->set_coeff(-tmp);
		gain *= 0.5*(1+tmp);
	  }
      for (j=odd;j<n2;j++) {                       
		td = a*a - 2*a*roots[j].real() + magsq(roots[j]);
		roots[j] = complex<double>((a*a - magsq(roots[j])),
					2.0*a*roots[j].imag())/td;
	  }
	}
	//! Get 2nd order IIR coefficients
	void get_coeff(long n2) {
        for (int j=odd;j<n2;j++) { 
		  iir[j-odd].set_coeff(-2*roots[j].real(),magsq(roots[j]));
		  gain *= (magsq(roots[j]) - 2*roots[j].real() + 1.0)/(4*magsq(roots[j])); 
	  }
	}
};
} // namespace SPUC 
