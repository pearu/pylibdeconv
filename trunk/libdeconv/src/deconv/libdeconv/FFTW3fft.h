/*
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * Author:    Yuansheng Sun (yuansheng-sun@uiowa.edu)
 * Copyright: University of Iowa 2006
 *
 * Filename:  FFTW3fft.h
 */


#ifndef FFTW3FFT_H
#define FFTW3FFT_H


#include <boost/shared_ptr.hpp>
#include <fftw3.h>


#define FFTW3_FLAG FFTW_MEASURE


/* 
	This class provides basic FFT routines developped based on FFTW3.
	It only supports 3-D complex-to-real (c2r, forward) or real-to-complex (r2c, backward) transforms 
	and is dedicated for 3-D deconvolution algorithms.

	Use it in 2 steps :
	- step 1 : create  a plan -> FFTW3_FFT::Ptr plan( new FFTW3_FFT( DimX, DimY, DimZ, IsForward, IsDouble, status )
	- step 2 : execute a plan -> p->execute( buf1, buf2, buf3 )

	<DimX> is the fastest varying dimension of a transform.
	<DimZ> are the slowest varying dimension of a transform.
	<_FFTsize> = DimX*DimY*(DimZ/2+1) if <DimZ> is even.
	<_FFTsize> = DimX*DimY*(DimZ+1)/2 if <DimZ> is odd.

	<IsDouble> indicates whether the floating type of the transformed data is double or single.

	<IsForward> = true : r2c transform plan, under this case
	<buf1> is input real, <buf2> is output real and <buf3> is output imaginary.
	<status> = 1 : <buf1> must be different than <buf2> and <buf3>, all arrays have the same size of DimX*DimY*DimZ. 
	<status> = 2 : real arrays overlap, buf1 and <buf2> must be same, all arrays have the same size of DimX*DimY*DimZ.
	<status> = 3 : <buf1> has the size of DimX*DimY*DimZ, <buf2> and <buf3> have the same size of _FFTsize.
  
	<IsForward> = false: c2r transform plan, under this case 
	<buf1> is input real and <buf2> is input imaginary, <buf3> is output real.
	<status> = 1 : <buf3> must be different than <buf1> and <buf2>, all arrays have the same size of DimX*DimY*DimZ 
	<status> = 2 : real arrays overlap, <buf3> and <buf1> must be same, all arrays have the same size of DimX*DimY*DimZ.
	<status> = 3 : <buf3> has the size of DimX*DimY*DimZ, <buf1> and <buf2> have the same size of <_FFTsize>.
	
	All date arrays must be one-dimensional and data is stored as "x + y*DimX + z*DimY*DimX".
	
	Throw: throw an error if fail.
*/
class FFTW3_FFT
{
	public:
	typedef boost::shared_ptr< FFTW3_FFT > Ptr ;

	~FFTW3_FFT() ;

	FFTW3_FFT( int DimX, int DimY, int DimZ, bool IsForward, bool IsDouble, int status ) ;

	void execute( double * buf1, double * buf2, double * buf3 ) ;
	void execute( float  * buf1, float  * buf2, float  * buf3 ) ;

	int   DimX()       { return _DimX ;      }
	int   DimY()       { return _DimY ;      }
	int   DimZ()       { return _DimZ ;      }
	int   FFTsize()    { return _FFTsize ;   }
	int   status()     { return _status ;    }
	bool  IsDouble()   { return _IsDouble ;  }
	bool  IsForward()  { return _IsForward ; }

        
	protected:
	int         _DimX ;
	int         _DimY ;
	int         _DimZ ;
	int         _FFTsize ;
	int         _status ;
	bool        _IsDouble ;
	bool        _IsForward ;
	double      _weight ;
	fftw_plan   _dplan ;
	fftwf_plan  _splan ;
} ;



/*
	This function provides a forward r2c FFT routine developped based on FFTW3 and 
	is dedicated for 3-D deconvolution algorithms. 
	Both single and double floating data types are supported. 
  	
	<DimX> is the fastest varying dimension of a transform.
	<DimZ> are the slowest varying dimension of a transform.
  	
	<in> is the input real.
	<out_re> and <out_im> are output real and imaginary, they must be different and not NULL.
	<in> and <out> can not be overlapped.  	
	All date arrays must be one-dimensional and data is stored as "x + y*DimX + z*DimY*DimX".
  	
	Throw: throw an error if fail.
*/
void fft3d( int DimX, int DimY, int DimZ, double * in, double * out_re, double * out_im ) ;
void fft3d( int DimX, int DimY, int DimZ, float  * in, float  * out_re, float  * out_im ) ;



/* 
	The following functions provide 1D/2D/3D FFT routines developped based on FFTW3.
	Both single and double floating data types are supported. 

	fft3d - 3D transform
	fft2d - 2D transform
	fft1d - 1D transform 
	
	<DimX> is the fastest varying dimension of a transform.
	<DimZ> are the slowest varying dimension of a transform.

	<IsForward> is true to indicate a forward transform, otherwise it is a backward transform. 

	<IsShift> is true to shift DC to the center in the frequency domain in a forward transform
	or you can also apply the same shift routine in a backward transform.

	<in_re>  and <in_im>  are input  real and imaginary, they must be different and not NULL.
	<out_re> and <out_im> are output real and imaginary, they must be different and not NULL.
	<in> and <out> can be overlapped. 
	All data arrays must be one-dimensional and data is stored as "x + y*DimX + z*DimY*DimX".
			    
	Throw: throw an error if fail.
*/
void fft3d( int DimX, int DimY, int DimZ, bool IsForward, bool IsShift, 
            double * in_re, double * in_im, double * out_re, double * out_im ) ;

void fft3d( int DimX, int DimY, int DimZ, bool IsForward, bool IsShift, 
            float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;

void fft2d( int DimX, int DimY, bool IsForward, bool IsShift, 
            double * in_re, double * in_im, double * out_re, double * out_im ) ; 

void fft2d( int DimX, int DimY, bool IsForward, bool IsShift, 
            float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;

void fft1d( int DimX, bool IsForward, bool IsShift,
            double * in_re, double * in_im, double * out_re, double * out_im ) ; 

void fft1d( int DimX, bool IsForward, bool IsShift,
            float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;


#endif  /*   #include "FFTW3fft.h"   */
