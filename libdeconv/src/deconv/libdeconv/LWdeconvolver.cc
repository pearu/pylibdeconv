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
 * Filename:  LWdeconvolver.cc
 */


#include <math.h>
#include "LWdeconvolver.h"

/* public functions */

LWdeconvolver::LWdeconvolver()
{
	init() ; 
}
	
LWdeconvolver::~LWdeconvolver()
{
	if (_FFTplanf)
		delete _FFTplanf ;
		
	if (_FFTplanb)
		delete _FFTplanb ;
}

void LWdeconvolver::init( bool IsApplyNorm, bool IsTrackLike, bool IsTrackMax, bool IsCheckStatus )
{
	_setControlFlags( IsCheckStatus, IsApplyNorm, IsTrackMax, IsTrackLike ) ;
	
	setConditioningIteration() ;
	setConditioningValue() ;
	setConditioningTolerance() ;
	
	setMaxRunIteration() ;
	setCriterion() ;
	
	_StartRunTime = 0 ;
	_StopRunTime  = 0 ;
	
	_DimX  = 0 ;
	_DimY  = 0 ;
	_DimZ  = 0 ;
	_Space = 0 ;
}

void LWdeconvolver::exportLW( const char * filename )
{
	FILE * fp = fopen( filename, "a+" ) ;
	if( fp )
	{
		_exportLWCG( fp ) ;
		fclose( fp ) ;
	}
	else
	{
		throw ErrnoError( std::string(filename) ) ;
	}
}

void LWdeconvolver::run( int DimX, int DimY, int DimZ, double * object_re, double * object_im, double * object, 
                         LWdws & ws, unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	double max_intensity = 0.0, cri = 1.0E+37, temp1, temp2 ;
	
	/* initialize running */
	_LWstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _LWprintStatus( 1 ) ;
	_initPSF( ws.size, object_im, ws.psf_re, ws.psf_im, FrequencySupport, ws.otf ) ;
	_initIMG( max_intensity, object_re, object, SpacialSupport ) ;
	fft3d( _DimX, _DimY, _DimZ, object_re, ws.image_re, ws.image_im ) ;	
	if( _CheckStatus ) _LWprintStatus( 2 ) ;
			
	/* start conditioning */
	if( _ConditioningIteration > 0 )
	{
		if( _CheckStatus ) _LWprintStatus( 3 ) ;
		for( int i = 0 ; i < _Space ; i++ ) ws.object0[i] = object[i] ;
		_ConditioningValue = _runConditioning( ws.size, ws.object0, object, object_re, object_im,
		                     ws.image_re, ws.image_im, ws.psf_re, ws.psf_im, ws.otf, SpacialSupport ) ;
		for( int i = 0 ; i < _Space ; i++ ) object[i] = ws.object0[i] ;
		delete [] ws.object0 ;
		ws.object0 = NULL ;
		if( _CheckStatus ) _LWprintStatus( 4 ) ;
	}
	
	/* initialize arrays in deconvolution loop */
	if( !_TrackLikelihood )
	{
		for( int i = 0 ; i < ws.size ; i++ )
		{
			         temp1 = ws.otf[i] + _ConditioningValue ;
			         temp2 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_im[i] = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_re[i] = temp2 ;
			     ws.otf[i] = ws.otf[i] / temp1 ;
		}
	}
	
	/* deconvolution loop */
	if( _CheckStatus ) _LWprintStatus( 5 ) ;
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _LWprintStatus( 6 ) ;		
		_FFTplanf->execute( object, object_re, object_im  ) ;		
		if( _TrackLikelihood )
		{
			_Likelihood.push_back( _getLikelihood( ws.size, ws.image_re, ws.image_im, 
			                       ws.psf_re, ws.psf_im, object_re, object_im ) ) ;
     			_LWupdate1( object_re, object_im, ws ) ;   
		}
		else
		{
			_LWupdate2( object_re, object_im, ws ) ; 
		}		
		for( int i = 0 ; i < _Space ; i++ )
		{
			object_im[i]  = object[i] ;
			   object[i] += object_re[i] ;
			if( object[i] < 0.0 ) object[i] = 0.0 ;
		}     		
		_getUpdate( object, object_im, SpacialSupport ) ;   
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}
		if( _CheckStatus ) _LWprintStatus( 7 ) ;
	}	
	
	/* end deconvolution */
	_LWfinishRun( ws ) ;	
}

void LWdeconvolver::run( int DimX, int DimY, int DimZ, float * object_re, float * object_im, float * object, 
                         LWsws & ws, unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	float  max_intensity = 0.0, temp1, temp2 ;
	double cri = 1.0E+37 ;

	/* initialize running */
	_LWstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _LWprintStatus( 1 ) ;
	_initPSF( ws.size, object_im, ws.psf_re, ws.psf_im, FrequencySupport, ws.otf ) ;
	_initIMG( max_intensity, object_re, object, SpacialSupport ) ;
	fft3d( _DimX, _DimY, _DimZ, object_re, ws.image_re, ws.image_im ) ;	
	if( _CheckStatus ) _LWprintStatus( 2 ) ;
	
	/* start conditioning */
	if( _ConditioningIteration > 0 )
	{
		if( _CheckStatus ) _LWprintStatus( 3 ) ;
		for( int i = 0 ; i < _Space ; i++ ) ws.object0[i] = object[i] ;
		_ConditioningValue = _runConditioning( ws.size, ws.object0, object, object_re, object_im, 
		                     ws.image_re, ws.image_im, ws.psf_re, ws.psf_im, ws.otf, SpacialSupport ) ;
		for( int i = 0 ; i < _Space ; i++ ) object[i] = ws.object0[i] ;
		delete [] ws.object0 ;
		ws.object0 = NULL ;
		if( _CheckStatus ) _LWprintStatus( 4 ) ;
	}
	
	/* initialize arrays in deconvolution loop */
	if( !_TrackLikelihood )
	{
		for( int i = 0 ; i < ws.size ; i++ )
		{
			         temp1 = ws.otf[i] + _ConditioningValue ;
			         temp2 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_im[i] = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_re[i] = temp2 ;
			     ws.otf[i] = ws.otf[i] / temp1 ;
		}
	}
	
	/* deconvolution loop */
	if( _CheckStatus ) _LWprintStatus( 5 ) ;
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _LWprintStatus( 6 ) ;		
		_FFTplanf->execute( object, object_re, object_im  ) ;		
		if( _TrackLikelihood )
		{
			_Likelihood.push_back( _getLikelihood( ws.size, ws.image_re, ws.image_im, 
			                       ws.psf_re, ws.psf_im, object_re, object_im ) ) ;
			_LWupdate1( object_re, object_im, ws ) ;  
		}
		else
		{
			_LWupdate2( object_re, object_im, ws ) ; 
		}		
		for( int i = 0 ; i < _Space ; i++ )
		{
			object_im[i]  = object[i] ;
			   object[i] += object_re[i] ;
			if( object[i] < 0.0 ) object[i] = 0.0 ;
		}
		_getUpdate( object, object_im, SpacialSupport ) ; 
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}
		if( _CheckStatus ) _LWprintStatus( 7 ) ;
	}
	
	/* end deconvolution */
	_LWfinishRun( ws ) ;	
}

/* private functions */

void LWdeconvolver::_LWprintStatus( int stage )
{
	switch( stage )
	{	
		case 1:
			time( &_t0 ) ;
			std::cout << " LWdeconvolver::run starts initialization ... \n" ;
			break ;

		case 2:
			time( &_t1 ) ;
			std::cout << " --> Calculate FFT on the input PSF.\n" ;
			if( _ApplyFrequencySupport )
			{
				std::cout << " --> apply frequency support on the FFT of the input PSF.\n" ;
			}
			if( _ApplyNormalization )
			{
				std::cout << " --> normalize the input image and the first estimated object.\n" ;
				std::cout << " --> Calculate FFT on the normalized input image.\n" ;
			}
			else
			{
				std::cout << " --> Calculate FFT on the input image.\n" ;
			}
			std::cout << " LWdeconvolver::run completes initialization, elapsed "
			          << difftime( _t1, _t0 ) << " seconds.\n" ;
			if( _ApplySpacialSupport )
			{
				std::cout << " LWdeconvolver::run applys spacial support on the deconvolved object iteratively.\n" ;
			}
			break ;

		case 3:
			time( &_t0 ) ;
			std::cout << " LWdeconvolver::run starts conditioning ... \n" ;
			break ;

		case 4:
			time( &_t1 ) ;
			std::cout << " LWdeconvolver::run completes conditioning, elapsed "
			          << difftime( _t1, _t0 ) << " seconds.\n" ;
			break ;

		case 5:
			std::cout << " LWdeconvolver::run starts the loop with conditioning value = " 
			          << _ConditioningValue << ".\n" ;
			break ;

		case 6:
			time( &_t0 ) ;
			break ;
		
		case 7:
			time ( &_t1 ) ;
			if( _TrackLikelihood )
			{
				printf( " --> Iteration %4d -> Update = %12.6e , Likelihood = %12.6e", (int)_Update.size(), 
				        _Update[ _Update.size()-1 ], _Likelihood[ _Likelihood.size()-1 ] ) ;
				std::cout << " , elapsed " << difftime( _t1, _t0 ) << " seconds.\n" ;
                        }
			else
			{
				printf( " --> Iteration %4d -> Update = %12.6e", 
				         (int)_Update.size(), _Update[ _Update.size()-1 ] ) ;
				std::cout << " , elapsed " << difftime( _t1, _t0 ) << " seconds.\n" ;
			}
			break ;

		default:
			break ;
	}
}

void LWdeconvolver::_LWstartRun( int DimX, int DimY, int DimZ, LWdws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double)_Space * 3.0 ;

	time( &_StartRunTime ) ;
	std::cout << " LWdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " LWdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " LWdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " LWdeconvolution stopping criterion : " << _Criterion << "\n" ;

	time( &_t0 ) ;
	std::cout << " LWdeconvolver::run starts creating FFT plans ... \n" ;

	_FFTplanf = new FFTW3_FFT (_DimX, _DimY, _DimZ, true,  true, 1) ; 	
	_FFTplanb = new FFTW3_FFT (_DimX, _DimY, _DimZ, false, true, 2); 

	time( &_t1 ) ;
	std::cout << " LWdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;

	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 5.0 ) ;
	
	ws.psf_re   = new double[ ws.size ] ;
	ws.psf_im   = new double[ ws.size ] ;
	ws.image_re = new double[ ws.size ] ;
	ws.image_im = new double[ ws.size ] ;
	ws.otf      = new double[ ws.size ] ;
	
	if( _ConditioningIteration > 0 )
	{
		ws.object0 = new double[ _Space ] ;
		memory += (double)_Space ;
	}
	else    ws.object0 = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;

	std::cout << " LWdeconvolver::run are using " << (memory/1024.0/128.0) << " Mbytes memory.\n" ;
}

void LWdeconvolver::_LWstartRun( int DimX, int DimY, int DimZ, LWsws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double)_Space * 3.0 ;

	time( &_StartRunTime ) ;
	std::cout << " LWdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " LWdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " LWdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " LWdeconvolution stopping criterion : " << _Criterion << "\n" ;

	time( &_t0 ) ;
	std::cout << " LWdeconvolver::run starts creating FFT plans ... \n" ;

	_FFTplanf = new FFTW3_FFT (_DimX, _DimY, _DimZ, true,  false, 1 ); 
	
	_FFTplanb = new FFTW3_FFT (_DimX, _DimY, _DimZ, false, false, 2) ; 
		
	time( &_t1 ) ;
	std::cout << " LWdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;

	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 5.0 ) ;

	ws.psf_re   = new float[ ws.size ] ;
	ws.psf_im   = new float[ ws.size ] ;
	ws.image_re = new float[ ws.size ] ;
	ws.image_im = new float[ ws.size ] ;
	ws.otf      = new float[ ws.size ] ;
	if( _ConditioningIteration > 0 )
	{
		ws.object0 = new float[ _Space ] ;
		memory += (double)_Space ;
	}
	else    ws.object0 = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;

	std::cout << " LWdeconvolver::run are using " << (memory/1024.0/256.0) << " Mbytes memory.\n" ;   	
}

void LWdeconvolver::_LWfinishRun( LWdws & ws )
{
	if( ws.psf_re   != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im   != NULL ) delete [] ws.psf_im ;
	if( ws.image_re != NULL ) delete [] ws.image_re ;
	if( ws.image_im != NULL ) delete [] ws.image_im ;
	if( ws.otf      != NULL ) delete [] ws.otf ;
		
	time( &_StopRunTime ) ;
	std::cout << " LWdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}  

void LWdeconvolver::_LWfinishRun( LWsws & ws )
{
	if( ws.psf_re   != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im   != NULL ) delete [] ws.psf_im ;
	if( ws.image_re != NULL ) delete [] ws.image_re ;
	if( ws.image_im != NULL ) delete [] ws.image_im ;
	if( ws.otf      != NULL ) delete [] ws.otf ;
	
	time( &_StopRunTime ) ;
	std::cout << " LWdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}

void LWdeconvolver::_LWupdate1( double * object_re, double * object_im, LWdws & ws  )
{
	double temp1, temp2, temp3, temp4 ;
	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = ws.otf[i] + _ConditioningValue ;
		       temp2 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
		       temp3 = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
		       temp4 = ws.otf[i] / temp1 ;
		object_re[i] = temp2 - temp4 * object_re[i] ;
		object_im[i] = temp3 - temp4 * object_im[i] ;
	}
	       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;
}



void LWdeconvolver::_LWupdate1( float * object_re, float * object_im, LWsws & ws )
{
	float temp1, temp2, temp3, temp4 ;
	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = ws.otf[i] + _ConditioningValue ;
		       temp2 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
		       temp3 = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
		       temp4 = ws.otf[i] / temp1 ;
		object_re[i] = temp2 - temp4 * object_re[i] ;
		object_im[i] = temp3 - temp4 * object_im[i] ;
	}
       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;	
}



void LWdeconvolver::_LWupdate2( double * object_re, double * object_im, LWdws & ws )
{
	for( int i = 0 ; i < ws.size ; i++ )
	{
		object_re[i] = ws.image_re[i] - ws.otf[i] * object_re[i] ;
		object_im[i] = ws.image_im[i] - ws.otf[i] * object_im[i] ;
	}
       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;
}


 
void LWdeconvolver::_LWupdate2( float * object_re, float * object_im, LWsws & ws )
{
	for( int i = 0 ; i < ws.size ; i++ )
	{
		object_re[i] = ws.image_re[i] - ws.otf[i] * object_re[i] ;
		object_im[i] = ws.image_im[i] - ws.otf[i] * object_im[i] ;
	}
       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;
}
