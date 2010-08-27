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
 * Filename:  FFTW3fft.cc
 */


#include "MYdef.h"
#include "MYerror.h"
#include "SHIFTfft.h"
#include "FFTW3fft.h"


void FFTW3_set_1D_dims( boost::shared_array< fftw_iodim >& dims, int DimX )
{
	dims[0].n  = DimX ;
	dims[0].is = 1 ;
	dims[0].os = 1 ;
}



void FFTW3_set_2D_dims( boost::shared_array< fftw_iodim >& dims, int DimX, int DimY )
{
	dims[1].n  = DimY ;
	dims[1].is = DimX ;
	dims[1].os = DimX ;
	dims[0].n  = DimX ;
	dims[0].is = 1 ;
	dims[0].os = 1 ;
}



void FFTW3_set_3D_dims( boost::shared_array< fftw_iodim >& dims, int DimX, int DimY, int DimZ )
{
	dims[2].n  = DimZ ;
	dims[2].is = DimX * DimY ;
	dims[2].os = DimX * DimY ;
	dims[1].n  = DimY ;
	dims[1].is = DimX ;
	dims[1].os = DimX ;
	dims[0].n  = DimX ;
	dims[0].is = 1 ;
	dims[0].os = 1 ;
}



class FFTW3Error : public Error
{
	public:
	FFTW3Error( int status, bool IsForward = true )
	{
		switch( status )
		{
			case -1:
				_error << " FFTW3-FFT error : FFT plan has not been created.\n" ;
				break ;
		
			case -2:
				if( IsForward )
				{
					_error << " FFTW3-FFT error : wrong use of FFTW3_guru_split_r2c single plan with double data.\n" ;
				}
				else
				{
					_error << " FFTW3-FFT error : wrong use of FFTW3_guru_split_c2r single plan with double data.\n" ;
				}
				break ;

			case -3:
				if( IsForward )
				{
					_error << " FFTW3-FFT error : wrong use of FFTW3_guru_split_r2c double plan with single data.\n" ;
				}
				else
				{
					_error << " FFTW3-FFT error : wrong use of FFTW3_guru_split_c2r double plan with single data.\n" ;
				}
				break ;

			case 0:
				_error << " FFTW3-FFT error : construction of FFTW3 plan failed.\n" ;
				break ;

			default:
				_error << " FFTW3-FFT error : unknow error.\n" ;
				break ;
		}
	}
} ;



FFTW3_FFT::~FFTW3_FFT()
{
	if( _status )
	{
		if( _IsDouble )
		{
			fftw_destroy_plan( _dplan ) ;
		}
		else
		{
			fftwf_destroy_plan( _splan ) ;
		}
	}
}



FFTW3_FFT::FFTW3_FFT( int DimX, int DimY, int DimZ, bool IsForward, bool IsDouble, int status )
{
	_DimX = DimX ;
	_DimY = DimY ;
	_DimZ = DimZ ;

	if( DimZ % 2 == 0 ) _FFTsize = (DimZ/2 + 1) * DimY * DimX ;
	else                _FFTsize = ((DimZ+1)/2) * DimY * DimX ;
	
	_weight = (double) (DimX * DimY * DimZ) ;
	_IsForward = IsForward ;
	_IsDouble = IsDouble ;
	_status = status ;
	
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[3] ) ;
	FFTW3_set_3D_dims( dims, DimX, DimY, DimZ ) ;

	if( status == 1 )
	{
		if( IsDouble )
		{
			DoubleArray buf1( new double[ DimX * DimY * DimZ ] ) ;
			DoubleArray buf2( new double[ DimX * DimY * DimZ ] ) ;
			DoubleArray buf3( new double[ DimX * DimY * DimZ ] ) ;
			
			if( IsForward )
			{
				_dplan = fftw_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
			else
			{	
				_dplan = fftw_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
		}
		else
		{
			SingleArray buf1( new float[ DimX * DimY * DimZ ] ) ;
			SingleArray buf2( new float[ DimX * DimY * DimZ ] ) ;
			SingleArray buf3( new float[ DimX * DimY * DimZ ] ) ;
			
			if( IsForward )
			{
				_splan = fftwf_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
			else
			{	
				_splan = fftwf_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
		}
	}
	else if( status == 2 )
	{	
		if( IsDouble )
		{
			DoubleArray buf1( new double[ DimX * DimY * DimZ ] ) ;
			DoubleArray buf2( new double[ DimX * DimY * DimZ ] ) ;
			
			if( IsForward )
			{
				_dplan = fftw_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf1.get(), buf2.get(), FFTW3_FLAG ) ;
			}
			else
			{	
				_dplan = fftw_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf1.get(), FFTW3_FLAG ) ;
			}
		}
		else
		{
			SingleArray buf1( new float[ DimX * DimY * DimZ ] ) ;
			SingleArray buf2( new float[ DimX * DimY * DimZ ] ) ;
			
			if( IsForward )
			{
				_splan = fftwf_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf1.get(), buf2.get(), FFTW3_FLAG ) ;
			}
			else
			{	
				_splan = fftwf_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf1.get(), FFTW3_FLAG ) ;
			}
		}
	}
	else if( status == 3 )
	{
		if( IsDouble )
		{
			if( IsForward )
			{
				DoubleArray buf1( new double[ DimX * DimY * DimZ ] ) ;
				DoubleArray buf2( new double[ _FFTsize ] ) ;
				DoubleArray buf3( new double[ _FFTsize ] ) ;
				
				_dplan = fftw_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
			else
			{	
				DoubleArray buf1( new double[ _FFTsize ] ) ;
				DoubleArray buf2( new double[ _FFTsize ] ) ;
				DoubleArray buf3( new double[ DimX * DimY * DimZ ] ) ;
				
				_dplan = fftw_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
		}
		else
		{
			if( IsForward )
			{
				SingleArray buf1( new float[ DimX * DimY * DimZ ] ) ;
				SingleArray buf2( new float[ _FFTsize ] ) ;
				SingleArray buf3( new float[ _FFTsize ] ) ;
				
				_splan = fftwf_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
			else
			{					
				SingleArray buf1( new float[ _FFTsize ] ) ;
				SingleArray buf2( new float[ _FFTsize ] ) ;
				SingleArray buf3( new float[ DimX * DimY * DimZ ] ) ;
				
				_splan = fftwf_plan_guru_split_dft_c2r( 3, dims.get(), 0, NULL, buf1.get(), buf2.get(), buf3.get(), FFTW3_FLAG ) ;
			}
		}
	}
	else 
	{
		throw FFTW3Error( 0 ) ;
	}
}


	
void FFTW3_FFT::execute( double * buf1, double * buf2, double * buf3 )
{
	if( _IsDouble )
	{
		if( _IsForward )
		{
			fftw_execute_split_dft_r2c( _dplan, buf1, buf2, buf3 ) ;
		}
		else
		{
			fftw_execute_split_dft_c2r( _dplan, buf1, buf2, buf3 ) ;

			for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ )
			{
				buf3[i] /= _weight ;
			}
		}
	}
	else
	{
		throw FFTW3Error( -2, _IsForward ) ;
	}
}



void FFTW3_FFT::execute( float * buf1, float * buf2, float * buf3 )
{
	if( !_IsDouble )
	{
		if( _IsForward )
		{
			fftwf_execute_split_dft_r2c( _splan, buf1, buf2, buf3 ) ;
		}
		else
		{
			fftwf_execute_split_dft_c2r( _splan, buf1, buf2, buf3 ) ;

			for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ )
			{
				buf3[i] /= _weight ;
			}
		}
	}
	else
	{
		throw FFTW3Error( -3, _IsForward ) ;
	}
}



void fft3d( int DimX, int DimY, int DimZ, double * in, double * out_re, double * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[3] ) ;
	FFTW3_set_3D_dims( dims, DimX, DimY, DimZ ) ;	
	
	fftw_plan p = fftw_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, in, out_re, out_im, FFTW_ESTIMATE ) ;

	if( p )
	{
		fftw_execute_split_dft_r2c( p, in, out_re, out_im ) ;
	}
	else
	{
		throw FFTW3Error( 0 ) ;		
	}
}



void fft3d( int DimX, int DimY, int DimZ, float * in, float * out_re, float * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[3] ) ;
	FFTW3_set_3D_dims( dims, DimX, DimY, DimZ ) ;	
	
	fftwf_plan p = fftwf_plan_guru_split_dft_r2c( 3, dims.get(), 0, NULL, in, out_re, out_im, FFTW_ESTIMATE ) ;

	if( p )
	{
		fftwf_execute_split_dft_r2c( p, in, out_re, out_im ) ;
	}
	else
	{
		throw FFTW3Error( 0 ) ;		
	}
}

		

void fft3d( int DimX, int DimY, int DimZ, bool IsForward, bool IsShift, 
            double * in_re, double * in_im, double * out_re, double * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[3] ) ;
	FFTW3_set_3D_dims( dims, DimX, DimY, DimZ ) ;	

	if( IsForward ) 
	{
		fftw_plan p = fftw_plan_guru_split_dft( 3, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift3d( DimX, DimY, DimZ, in_re, in_im, in_re, in_im ) ;
			}

			fftw_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift3d( DimX, DimY, DimZ, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftw_plan p = fftw_plan_guru_split_dft( 3, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) ( DimX * DimY * DimZ ) ;

			fftw_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			double * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift3d( DimX, DimY, DimZ, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX * DimY * DimZ ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}



void fft3d( int DimX, int DimY, int DimZ, bool IsForward, bool IsShift, 
            float * in_re, float * in_im, float * out_re, float * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[3] ) ;
	FFTW3_set_3D_dims( dims, DimX, DimY, DimZ ) ;

	if( IsForward ) 
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 3, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift3d( DimX, DimY, DimZ, in_re, in_im, in_re, in_im ) ;
			}

			fftwf_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift3d( DimX, DimY, DimZ, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 3, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) ( DimX * DimY * DimZ ) ;

			fftwf_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			float * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift3d( DimX, DimY, DimZ, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX * DimY * DimZ ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}



void fft2d( int DimX, int DimY, bool IsForward, bool IsShift, 
            double * in_re, double * in_im, double * out_re, double * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[2] ) ;
	FFTW3_set_2D_dims( dims, DimX, DimY ) ;

	if( IsForward ) 
	{
		fftw_plan p = fftw_plan_guru_split_dft( 2, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift2d( DimX, DimY, in_re, in_im, in_re, in_im ) ;
			}

			fftw_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift2d( DimX, DimY, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftw_plan p = fftw_plan_guru_split_dft( 2, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) ( DimX * DimY ) ;

			fftw_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			double * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift2d( DimX, DimY, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX * DimY ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}



void fft2d( int DimX, int DimY, bool IsForward, bool IsShift, 
            float * in_re, float * in_im, float * out_re, float * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[2] ) ;
	FFTW3_set_2D_dims( dims, DimX, DimY ) ;	

	if( IsForward ) 
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 2, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift2d( DimX, DimY, in_re, in_im, in_re, in_im ) ;
			}

			fftwf_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift2d( DimX, DimY, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 2, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) ( DimX * DimY ) ;

			fftwf_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			float * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift2d( DimX, DimY, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX * DimY ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}



void fft1d( int DimX, bool IsForward, bool IsShift, 
            double * in_re, double * in_im, double * out_re, double * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[1] ) ;
	FFTW3_set_1D_dims( dims, DimX ) ;

	if( IsForward ) 
	{
		fftw_plan p = fftw_plan_guru_split_dft( 1, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift1d( DimX, in_re, in_im, in_re, in_im ) ;
			}

			fftw_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift1d( DimX, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftw_plan p = fftw_plan_guru_split_dft( 1, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) DimX ;

			fftw_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			double * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift1d( DimX, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}



void fft1d( int DimX, bool IsForward, bool IsShift, 
            float * in_re, float * in_im, float * out_re, float * out_im )
{
	boost::shared_array< fftw_iodim > dims( new fftw_iodim[1] ) ;
	FFTW3_set_1D_dims( dims, DimX ) ;

	if( IsForward ) 
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 1, dims.get(), 0, NULL, in_re, in_im, out_re, out_im, FFTW_ESTIMATE ) ;

		if( p )
		{
			if( IsShift )
			{
				shift1d( DimX, in_re, in_im, in_re, in_im ) ;
			}

			fftwf_execute_split_dft( p, in_re, in_im, out_re, out_im ) ;

			if( IsShift && in_re != out_re && in_im != out_im )
			{
				shift1d( DimX, in_re, in_im, in_re, in_im ) ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
	else
	{
		fftwf_plan p = fftwf_plan_guru_split_dft( 1, dims.get(), 0, NULL, in_im, in_re, out_im, out_re, FFTW_ESTIMATE ) ;

		if( p )
		{
			double weight = (double) DimX ;

			fftwf_execute_split_dft( p, in_im, in_re, out_im, out_re ) ;
			float * temp = out_im ;
			out_im = out_re ;
			out_re = temp ;

			if( IsShift )
			{
				shift1d( DimX, out_re, out_im, out_re, out_im ) ;
			}

			for( int i = 0 ; i < DimX ; i++ )
			{
				out_re[i] /= weight ;
				out_im[i] /= weight ;
			}
		}
		else
		{
			throw FFTW3Error( 0 ) ;		
		}
	}
}
