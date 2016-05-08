#if !defined( _CFFT_H_INC__ )
#define _CFFT_H_INC__ 1
#include <math.h>  // for sin and cos
namespace SPUC {
//!  \brief  General Purpose C++ complex FFT Transform
//!  \ingroup miscfunc
//
//! This is a general-purpose C++ complex FFT transform class.
//! it is defined as a template over a complex type. For instance,
//! if using gnu gcc, the complex type is
//!			complex<double>    
//! And you declare the cfft class as
//!          cfft<complex<double>>
//!
//! The underlying CPLX type requires:
//!  CPLX()
//!  operator = , CPLX(CPLX const&)
//!	CPLX(double,double)	[used on cos/sin]
//!	CPLX operator*( CPLX , double )
//!  CPLX conj(CPLX const &);		[conjugate]
//!  ComPlex::operator @ (CPLX , CPLX )	[ where @ =* + - ]
template <class CPLX> class cfft {
    int N, log2N;		//! these define size of FFT buffer
    CPLX *w;			//! array [N/2] of cos/sin values
    int *bitrev;		//! bit-reversing table, in 0..N
    double fscales[2];	//! f-transform scalings
    double iscales[2];	//! i-transform scales
    void fft_func( CPLX *buf, int iflag );
	
public:
    cfft( int size,		// size is power of 2
		  double scalef1 = 0.5, double scalef2 = 1.0,	// fwd transform scalings
		  double scalei1 = 1.0, double scalei2 = 1.0	// rev xform
		  );
    ~cfft();
    inline void fft(CPLX *buf)		//! perform forward fft on buffer
		{
			fft_func( buf, 0 );
		}
    inline void ifft(CPLX *buf)		//! perform reverse fft on buffer
		{
			fft_func( buf, 1 );
		}
    inline int length() const { return N; }
    //! used to fill in last half of complex spectrum of real signal
    //! when the first half is already there.
    void hermitian( CPLX *buf );

}; // class cfft
//////////////////////////////  cfft methods //////////////////////////////
/*
 * constructor takes an int, power-of-2.
 * scalef1,scalef2, are the post-pass and post-transform
 * scalings for the forward transform; scalei1 and scalei2 are
 * the same for the inverse transform.
 */
template <class CPLX>
cfft<CPLX>::cfft( int size, double scalef1, double scalef2,
		  double scalei1, double scalei2 )
{
    int i,j,k;
    double t;

    fscales[0] = scalef1;
    fscales[1] = scalef2;
    iscales[0] = scalei1;
    iscales[1] = scalei2;

    for( k = 0; ; ++k ){
	if( (1<<k) == size) break;
	if( k==18 || (1<<k) >size )
	    throw "cfft: size not power of 2";
    }
    N = 1<<k;
    log2N = k;

    bitrev = new int [N];

    if( k > 0 )
	w = new CPLX[ N>>1 ];
    else
	w = NULL;
    if( bitrev == NULL || ((k>0) && w == NULL) )
	throw "cfft: out of memory";
	
    //
    // do bit-rev table
    //
    bitrev[0] = 0;

    for( j = 1; j < N; j<<=1 ){
	for( i = 0; i < j ; ++i ){
	    bitrev[i] <<= 1;
	    bitrev[i+j] = bitrev[i]+1;
	}
    }
    //
    // prepare the cos/sin table. This is bit-reversed, and goes
    // like this: 0, 90, 45, 135, 22.5 ...  for N/2 entries.
    //
    if( k > 0 ){
	CPLX ww;
	k = (1<<(k-1));
	for( i = 0; i < k; ++i ){
	    t = double(bitrev[i<<1]) * PI /double(k);
	    ww = CPLX( cos(t), sin(t));
	    w[i] = conj(ww);		// force limiting of imag part if applic.
	    //cout << w[i] << "\n";
	}
    }
}


/*
 * destructor frees the memory
 */
template <class CPLX>
cfft<CPLX>::~cfft()
{
    delete [] bitrev;
    if( w!= NULL )
	delete [] w;
}


/*
 * hermitian() assumes the array has been filled in with values
 * up to the center and including the center point. It reflects these
 * conjugate-wise into the last half.
 */
template <class CPLX>
void
cfft<CPLX>::hermitian(CPLX *buf)
{
    int i,j;
    if( N <= 2 ) return; 	// nothing to do
    i = (N>>1)-1;			// input
    j = i+2;				// output
    while( i > 0 ){
	buf[j] = conj(buf[i]);
	--i;
	++j;
    }
}

	
/*
 * cfft::fft_func(buf,0) performs a forward fft on the data in the buffer specified.
 * cfft::fft_func(buf,1) performs an inverse fft on the data in the buffer specified.
 */
template <class CPLX>
void
cfft<CPLX>::fft_func( CPLX *buf, int iflag )
{
    int i,j,k;
    CPLX *buf0,*buf2,*bufe;
    CPLX z1,z2,zw;
    double *sp,s;

    sp = iflag ? iscales : fscales;
    s = sp[0];		// per-pass scale

    if( log2N == 0 ){		// only 1 element !
	s = sp[1];
	buf[0] = buf[0] * s;	// final scale only
	return;
    }
    //
    // first pass:
    //  1st element  = sum of 1st & middle, middle element = diff.
    // repeat N/2 times.

    k = N>>1;

    if( log2N == 1 )
	s *= sp[1];	// final scale
	
    buf2 = buf + k;
    for( i = 0; i < k; ++i ){		// first pass is faster
	z1 = buf[i] + buf2[i];
	z2 = buf[i] - buf2[i];
	buf[i] = z1 * s;
	buf2[i] = z2 * s;
    }
    if( log2N == 1 ) return;	// only 2!

    k>>=1;				// k is N/4 now
    bufe = buf+N;		// past end
    for( ; k; k>>=1 ){
	if( k == 1 ){	// last pass - include final scale 
	    s *= sp[1];	// final scale
	}
	buf0=buf;
	for( j = 0; buf0<bufe; ++j ){
	    if( iflag ){
		zw = conj(w[j]);
	    }else{
		zw = w[j];		/* get w-factor */
	    }
	    buf2 = buf0+k;
	    for( i = 0; i < k; ++i ){	// a butterfly
		z1 = zw * buf2[i];
		z2 = buf0[i] + z1;
		buf2[i] = (buf0[i] - z1)*s;
		buf0[i] = z2 * s;
	    }
	    buf0 += (k<<1);
	}
    }
    // bitrev the sucker 
    for( i = 0; i < N; ++i ){
	j = bitrev[i];
	if( i <= j ) continue;		// don't do these
	z1 = buf[i];
	buf[i] = buf[j];
	buf[j] = z1;
    }
}
////////////////////////////// end cfft //////////////////////////////
}
#endif
