/* *************************************************************************
ifft - In-place radix 2 decimation in time inverse FFT

Requires pointer to complex array, x and power of 2 size of FFT, m
(size of FFT = 2**m).  Places inverse FFT output on top of input
frequency domain COMPLEX array.

complex<double> ifft(complex<double> *x, int m)

*************************************************************************/
#include <complex.h>
#include <iostream>
#include <fstream>

#include "idft.h"

using namespace std;

namespace SPUC {
void ifft(complex<double> *y, int m)
{
  complex<double> u,temp,tm;
  int wptr;
  int i,j,k,l,le,windex;
  float scale;

  /*  n = 2**m = inverse fft length */
  int n = 1 << m;
  le = n/2;
  complex<double> *w = new complex<double> [n+1];
  complex<double> *x  = new complex<double> [n+1];

  /*  calculate the w values recursively */
  complex<double> w_inc = expj(TWOPI/n);
  complex<double> w_x = w_inc;

  for (j=0; j<le; j++) {
    w[j] = w_x;
    w_x *= w_inc;
	x[j] = y[j];
	x[j+le] = y[j+le];
  }

  /*  start inverse fft */
  le = n;
  windex = 1;
  for (l=0; l<m; l++) {
    le = le/2;
    /*  first iteration with no multiplies */
    for (i=0; i<n; i += 2*le) {
      temp = (x[i] + x[i+le]);
      x[i+le] = (x[i] - x[i+le]);
      x[i] = temp;
    }

    /*  remaining iterations use stored w */
    wptr = windex-1;
    for (j=1; j<le ; j++) {
      for (i=j; i<n; i += 2*le) {
		temp = x[i] + x[i+le];
		tm = x[i] - x[i+le];
		x[i+le] = tm*w[wptr];
		x[i] = temp;
      }
      wptr += windex;
    }
    windex = 2*windex;
  }            

  /*  rearrange data by bit reversing */

  j=0;
  for (i=1; i<(n-1); i++) {
    k = n/2;
    while(k <= j) {
      j = j - k;
      k = k/2;
    }
    j = j + k;
    if (i < j) {
      temp = x[j];
      x[j] = x[i];
      x[i] = temp;
    }
  }
  
  /*  scale all results by 1/n */
  scale = (float)(1.0/n);
  for (i=0; i<n; i++) y[i] = scale*x[i];

  delete [] x;
  delete [] w;

}
}
