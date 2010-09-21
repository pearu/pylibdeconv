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
 * Filename:  EMdeconvolver.cc
 */


#include <math.h>
#include "EMdeconvolver.h"

EMdeconvolver::EMdeconvolver() :deconvolver()
{ 
	init() ; 
}		

EMdeconvolver::~EMdeconvolver()
{
	if (_FFTplanf)
		delete _FFTplanf ;
		
	if (_FFTplanb)
		delete _FFTplanb ;
}

/* public functions */

void EMdeconvolver::setEMIRpenalty( double penalty )      
{ 
	if( penalty > EMDepsilon )
	{ 
		_EMIRpenalty = penalty ; 
	}
	else
	{
		_EMIRpenalty = -1.0 ;
		std::cout << " Intensity_Penalty was not properly set and will be calculated in EMdeconvolver::run().\n" ;
	}
}


	
void EMdeconvolver::init( bool IsAccelerate, bool IsApplyNorm, bool IsTrackLike, bool IsTrackMax, bool IsCheckStatus )
{
	_setControlFlags( IsCheckStatus, IsApplyNorm, IsTrackMax, IsTrackLike ) ;
	
	_Accelerate = IsAccelerate ;
	
	setEMIRiteration() ;
	_EMIRpenalty = -1.0 ;
	
	setMaxRunIteration() ;
	setCriterion() ;
	
	_StartRunTime = 0 ;
	_StopRunTime  = 0 ;
	
	_DimX  = 0 ;
	_DimY  = 0 ;
	_DimZ  = 0 ;
	_Space = 0 ;
}



void EMdeconvolver::exportEM( const char * filename )
{
	FILE * fp = fopen( filename, "a+" ) ;
	if( fp )
	{
		_exportCommon( fp ) ;
			
		if( _EMIRiteration > 0 )
		{
			fprintf( fp, "%d -> Apply intensity regularization every %d iterations in the deconvolution loop.\n", _EMIRiteration, _EMIRiteration ) ; 
			fprintf( fp, "\n" ) ;
		}
			
		if( _EMIRpenalty > 0.0 )
		{
			fprintf( fp, "%e -> Penalty value applied for intensity regularization in the deconvolution loop.\n", _EMIRpenalty ) ;
			fprintf( fp, "\n" ) ;
		}
			
		fprintf( fp, "%d -> Apply newton acceleration iteratively in the deconvolution loop.\n", ((int) _Accelerate) ) ; 
		fprintf( fp, "\n" ) ;
			
		fclose( fp ) ;
	}
	else
	{
		throw ErrnoError( std::string(filename) ) ;
	}
}



void EMdeconvolver::run( int DimX, int DimY, int DimZ, double * image, double * rat, double * object, EMdws & ws,
                         unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	double max_intensity = 0.0, cri = 1.0E+37 ;

	/* initialize running */
	_EMstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _EMprintStatus( 1 ) ;
	_initPSF( ws.size, rat, ws.psf_re, ws.psf_im, FrequencySupport ) ;
	_initIMG( max_intensity, image, object, SpacialSupport ) ;	
	if( _CheckStatus ) _EMprintStatus( 2 ) ;	
	
	/* start regularization */
	if( SpacialSupport != NULL ) _ApplySpacialSupport = true ;
	if( _EMIRiteration > 0 )
	{
		if( _EMIRpenalty < EMDepsilon )
		{
			if( _ApplyNormalization ) _EMIRpenalty = rat[0] ;
			else
			{
				max_intensity = image[0] ;
				for( int i = 0 ; i < _Space ; i++ )
				{
					if( image[i] > max_intensity ) max_intensity = image[i] ;
				}
				_EMIRpenalty = rat[0] / max_intensity ; 
			}
		}
		if( _CheckStatus ) _EMprintStatus( 3 ) ;
	}	

	/* deconvolution loop */
	if( _CheckStatus ) _EMprintStatus( 4 ) ;	
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _EMprintStatus( 5 ) ;
		if ( _Accelerate )
		{ 
			_EMupdate2( image, rat, object, ws ) ;
		}
		else
		{
			_EMupdate1( image, rat, object, ws ) ;        		
		} 
        	if( _EMIRiteration > 0 )
		{
			if( (_Update.size() + 1) % _EMIRiteration == 0 )
			{
				for( int i = 0 ; i < _Space ; i++ ) 
				{
					object[i] = ( -1.0 + sqrt(1.0 + 2.0 * _EMIRpenalty * object[i]) ) / _EMIRpenalty ;
				}
			}
		}
		_getUpdate( object, ws.buf, SpacialSupport ) ;
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}		
		if( _CheckStatus ) _EMprintStatus( 6 ) ;
	}	
	
	/* end deconvolution */
	_EMfinishRun( ws ) ;	
}



void EMdeconvolver::run( int DimX, int DimY, int DimZ, float * image, float * rat, float * object, EMsws & ws,
                         unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	float  max_intensity = 0.0 ;
	double cri = 1.0E+37 ;

	/* initialize running */
	_EMstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _EMprintStatus( 1 ) ;
	_initPSF( ws.size, rat, ws.psf_re, ws.psf_im, FrequencySupport ) ;
	_initIMG( max_intensity, image, object, SpacialSupport ) ;	
	if( _CheckStatus ) _EMprintStatus( 2 ) ;	
	
	/* start regularization */
	if( _EMIRiteration > 0 )
	{
		if( _EMIRpenalty < EMDepsilon )
		{
			if( _ApplyNormalization ) _EMIRpenalty = rat[0] ;
			else
			{
				max_intensity = image[0] ;
				for( int i = 0 ; i < _Space ; i++ )
				{
					if( image[i] > max_intensity ) max_intensity = image[i] ;
				}
				_EMIRpenalty = rat[0] / max_intensity ; 
			}
		}
		if( _CheckStatus ) _EMprintStatus( 3 ) ;
	}	

	/* deconvolution loop */
	if( _CheckStatus ) _EMprintStatus( 4 ) ;	
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _EMprintStatus( 5 ) ;
		if ( _Accelerate )
		{ 
			_EMupdate2( image, rat, object, ws ) ;
		}
		else
		{
			_EMupdate1( image, rat, object, ws ) ;        		
		} 
		if( _EMIRiteration > 0 )
		{
			if( (_Update.size() + 1) % _EMIRiteration == 0 )
			{
				for( int i = 0 ; i < _Space ; i++ ) 
				{
					object[i] = ( -1.0 + sqrt(1.0 + 2.0 * _EMIRpenalty * object[i]) ) / _EMIRpenalty ;
				}
			}
		}
		_getUpdate( object, ws.buf, SpacialSupport ) ;
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}		
		if( _CheckStatus ) _EMprintStatus( 6 ) ;
	}	
	
	/* end deconvolution */
	_EMfinishRun( ws ) ;	
}



/* private functions */

void EMdeconvolver::_EMprintStatus( int stage )
{
	switch( stage )
	{	
		case 1:
			time( &_t0 ) ;
			std::cout << " EMdeconvolver::run starts initialization ... \n" ;
			break ;

		case 2:
			time( &_t1 ) ;
			std::cout << " --> calculate FFT on the input PSF.\n" ;
			if( _ApplyFrequencySupport )
			{
				std::cout << " --> apply frequency support on the FFT of the input PSF.\n" ;
			}
			if( _ApplyNormalization )
			{
				std::cout << " --> normalize the input image and the first estimated object.\n" ;
			}                  
			std::cout << " EMdeconvolver::run completes initialization, elapsed "
			          << difftime( _t1, _t0 ) << " seconds.\n" ;
			if( _ApplySpacialSupport )
			{
				std::cout << " EMdeconvolver::run applys spacial support on the deconvolved object iteratively.\n" ;
			}
			break ;

		case 3:
			std::cout << " EMdeconvolver::run applys intensity regularization with penalty = " 
			          << _EMIRpenalty << " every " << _EMIRiteration << " iterations.\n" ;
			break ;

		case 4:
			std::cout << " EMdeconvolver::run starts the loop ...\n" ;
			break ;

		case 5:
			time( &_t0 ) ;
			break ;                	
		
		case 6:
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



void EMdeconvolver::_EMprintAcceleration( double alpha, double likelihood, int negative_points ) 
{
	printf( " --> Newton Acceleration Value = %9.6f -> Likelihood = %12.6e , Negative Points = %9d\n", 
	         alpha, likelihood, negative_points ) ; 
}



void EMdeconvolver::_EMprintAcceleration( double alpha ) 
{
	printf( " --> Newton Acceleration Value = %9.6f\n", alpha ) ; 
}



void EMdeconvolver::_EMstartRun( int DimX, int DimY, int DimZ, EMdws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double)_Space * 3 ;
		
	time( &_StartRunTime ) ;
	std::cout << " EMdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " EMdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " EMdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " EMdeconvolution stopping criterion : " << _Criterion << "\n" ;
        	
	time( &_t0 ) ;
	std::cout << " EMdeconvolver::run starts creating FFT plans ... \n" ;
        	
	_FFTplanf = new FFTW3_FFT (_DimX, _DimY, _DimZ, true,  true, 3) ;
	_FFTplanb = new FFTW3_FFT (_DimX, _DimY, _DimZ, false, true, 3) ;
		
	time( &_t1 ) ;
	std::cout << " EMdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;
	
	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 4.0 + (double)_Space ) ;
	
	ws.psf_re = new double[ ws.size ] ;
	ws.psf_im = new double[ ws.size ] ;
	ws.buf_re = new double[ ws.size ] ;
	ws.buf_im = new double[ ws.size ] ;
	ws.buf    = new double[ _Space  ] ;
	if( _Accelerate || _TrackLikelihood )
	{
		ws.eimg = new double[ _Space ] ;
		memory += ( (double)_Space ) ;
	}
	else    ws.eimg = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;
	
	std::cout << " EMdeconvolver::run are using " << (memory/1024.0/128.0) << " Mbytes memory.\n" ;
}



void EMdeconvolver::_EMstartRun( int DimX, int DimY, int DimZ, EMsws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double) _Space * 3.0 ;
		
	time( &_StartRunTime ) ;
	std::cout << " EMdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " EMdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " EMdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " EMdeconvolution stopping criterion : " << _Criterion << "\n" ;
        	
	time( &_t0 ) ;
	std::cout << " EMdeconvolver::run starts creating FFT plans ... \n" ;
        	
	_FFTplanf = new FFTW3_FFT (_DimX, _DimY, _DimZ, true,  false, 3) ;
	_FFTplanb = new FFTW3_FFT( _DimX, _DimY, _DimZ, false, false, 3) ;
		
	time( &_t1 ) ;
	std::cout << " EMdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;
	
	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 4.0 + (double)_Space ) ;
	
	ws.psf_re = new float[ ws.size ] ;
	ws.psf_im = new float[ ws.size ] ;
	ws.buf_re = new float[ ws.size ] ;
	ws.buf_im = new float[ ws.size ] ;
	ws.buf    = new float[ _Space  ] ;
	
	if( _Accelerate || _TrackLikelihood )
	{
		ws.eimg = new float[ _Space ] ;
		memory += ( (double)_Space ) ;
	}
	else    ws.eimg = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;
	
	std::cout << " EMdeconvolver::run are using " << (memory/1024.0/256.0) << " Mbytes memory.\n" ;
}



void EMdeconvolver::_EMfinishRun( EMdws & ws )
{
	if( ws.psf_re != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im != NULL ) delete [] ws.psf_im ;
	if( ws.buf_re != NULL ) delete [] ws.buf_re ;
	if( ws.buf_im != NULL ) delete [] ws.buf_im ;
	if( ws.buf    != NULL ) delete [] ws.buf ;
	if( ws.eimg   != NULL ) delete [] ws.eimg ;
		
	time( &_StopRunTime ) ;
	std::cout << " EMdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}  



void EMdeconvolver::_EMfinishRun( EMsws & ws )
{
	if( ws.psf_re != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im != NULL ) delete [] ws.psf_im ;
	if( ws.buf_re != NULL ) delete [] ws.buf_re ;
	if( ws.buf_im != NULL ) delete [] ws.buf_im ;
	if( ws.buf    != NULL ) delete [] ws.buf ;
	if( ws.eimg   != NULL ) delete [] ws.eimg ;
		
	time( &_StopRunTime ) ;
	std::cout << " EMdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}  



void EMdeconvolver::_EMupdate1( double * image, double * rat, double * object, EMdws & ws )
{
	double temp ;
	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;

	if( _TrackLikelihood )
	{
		double likelihood = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if ( rat[i] < EMDepsilon ) rat[i] = EMDepsilon ;
			likelihood += ( image[i] * log(rat[i]) - rat[i] ) ;
			rat[i] = image[i] / rat[i] ;
		}
		_Likelihood.push_back( likelihood ) ;
	}
	else
	{
		for( int i = 0 ; i < _Space ; i++ )
		{
			if ( rat[i] < EMDepsilon ) rat[i] = EMDepsilon ;
			rat[i] = image[i] / rat[i] ;
		}
	}
		
	_FFTplanf->execute( rat, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] + ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_im[i] * ws.psf_re[i] - ws.buf_re[i] * ws.psf_im[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;
		
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		ws.buf[i]  = object[i] ;	
		object[i] *= rat[i] ;
		if ( object[i] < 0.0 ) object[i] = 0.0 ;
	}
}



void EMdeconvolver::_EMupdate1( float * image, float * rat, float * object, EMsws & ws )
{
	float temp ;
	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;

	if( _TrackLikelihood )
	{
		double likelihood = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if ( rat[i] < EMSepsilon ) rat[i] = EMSepsilon ;
			likelihood += ( image[i] * log(rat[i]) - rat[i] ) ;
			rat[i] = image[i] / rat[i] ;
		}
		_Likelihood.push_back( likelihood ) ;
	}
	else
	{
		for( int i = 0 ; i < _Space ; i++ )
		{
			if ( rat[i] < EMSepsilon ) rat[i] = EMSepsilon ;
			rat[i] = image[i] / rat[i] ;
		}
	}
		
	_FFTplanf->execute( rat, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] + ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_im[i] * ws.psf_re[i] - ws.buf_re[i] * ws.psf_im[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;
		
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		ws.buf[i]  = object[i] ;	
		object[i] *= rat[i] ;
		if ( object[i] < 0.0 ) object[i] = 0.0 ;
	}
}



void EMdeconvolver::_EMupdate2( double * image, double * rat, double * object, EMdws & ws )
{
	double alpha = 1.0, alpha_new = 1.0, likelihood, temp, temp1, temp2 ;
	bool   continue_loop_1 = true, continue_loop_2 = true ;
	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, ws.eimg ) ;
     	
	for( int i = 0 ; i < _Space ; i++ )
	{
		if ( ws.eimg[i] < EMDepsilon ) ws.eimg[i] = EMDepsilon ;
		rat[i] = image[i] / ws.eimg[i] ;
	}
     	
	_FFTplanf->execute( rat, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] + ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_im[i] * ws.psf_re[i] - ws.buf_re[i] * ws.psf_im[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;
     	
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		ws.buf[i]  = object[i] ;
		object[i] *= ( rat[i] - 1.0 ) ;
	}
     	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;

	do
	{ 
		likelihood = 0.0 ;
		int j = 0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] < 0.0 && ( alpha_new * object[i] + ws.buf[i] ) < 0.0 ) j++ ;
			temp = ws.eimg[i] + alpha_new * rat[i] ;
			if ( temp > EMDepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
		}           			
		if( _CheckStatus ) _EMprintAcceleration( alpha_new, likelihood, j ) ;
		if( j > 0 )
		{
			while( alpha < alpha_new && continue_loop_1 )
			{
				alpha *= 1.5 ;
				int k = 0 ;
				likelihood = 0.0 ; 
				for( int i = 0 ; i < _Space ; i++ )
				{
					if( object[i] < 0.0 && ( alpha * object[i] + ws.buf[i] ) < 0.0 ) k++ ;
					temp = ws.eimg[i] + alpha * rat[i] ;
					if ( temp > EMDepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
				}
				if( _CheckStatus ) _EMprintAcceleration( alpha, likelihood, k ) ;       		
				if( k > 0 )
				{
					alpha /= 1.5 ;
					continue_loop_1 = false ;
				}
			}
			continue_loop_2 = false ;
		}
		if( continue_loop_2 )
		{
			alpha = alpha_new ;
			temp1 = 0.0 ;
			temp2 = 0.0 ;
			for( int i = 0 ; i < _Space ; i++ ) 
			{
				temp = ws.eimg[i] + alpha * rat[i] ; 
				if( fabs( temp ) > EMDepsilon )
				{ 
					temp1 += ( image[i] * rat[i] / temp - rat[i] ) ;
					temp2 += ( image[i] * rat[i] * rat[i] / ( temp * temp ) ) ;
				}
			}
			alpha_new = alpha + temp1 / temp2 ;
		}
	}
	while( fabs( alpha_new - alpha ) > alpha * 0.1 && continue_loop_2 ) ;
	
	if( _CheckStatus ) _EMprintAcceleration( alpha ) ;
             
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		object[i] = ws.buf[i] + alpha * object[i] ;
		if ( object[i] < 0.0 ) object[i] = 0.0 ;
	}
        
	if( _TrackLikelihood )
	{
		likelihood = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			temp = ws.eimg[i] + alpha * rat[i] ;
			if ( temp > EMDepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
		}
		_Likelihood.push_back( likelihood ) ;
	}
}



void EMdeconvolver::_EMupdate2( float * image, float * rat, float * object, EMsws & ws )
{
	float  temp ;
	double alpha = 1.0, alpha_new = 1.0, likelihood, temp1, temp2 ;
	bool   continue_loop_1 = true, continue_loop_2 = true ;
	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, ws.eimg ) ;
     	
	for( int i = 0 ; i < _Space ; i++ )
	{
		if ( ws.eimg[i] < EMSepsilon ) ws.eimg[i] = EMSepsilon ;
		rat[i] = image[i] / ws.eimg[i] ;
	}
     	
	_FFTplanf->execute( rat, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] + ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_im[i] * ws.psf_re[i] - ws.buf_re[i] * ws.psf_im[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;
     	
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		ws.buf[i]  = object[i] ;
		object[i] *= ( rat[i] - 1.0 ) ;
	}
     	
	_FFTplanf->execute( object, ws.buf_re, ws.buf_im ) ;
	for( int i = 0 ; i < ws.size ; i++ )
	{
		        temp = ws.buf_re[i] * ws.psf_re[i] - ws.buf_im[i] * ws.psf_im[i] ;  
		ws.buf_im[i] = ws.buf_re[i] * ws.psf_im[i] + ws.buf_im[i] * ws.psf_re[i] ;
		ws.buf_re[i] = temp ;
	}
	_FFTplanb->execute( ws.buf_re, ws.buf_im, rat ) ;

	do
	{ 
		likelihood = 0.0 ;
		int j = 0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] < 0.0 && ( alpha_new * object[i] + ws.buf[i] ) < 0.0 ) j++ ;
			temp = ws.eimg[i] + alpha_new * rat[i] ;
			if ( temp > EMSepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
		}           			
		if( _CheckStatus ) _EMprintAcceleration( alpha_new, likelihood, j ) ;
		if( j > 0 )
		{
			while( alpha < alpha_new && continue_loop_1 )
			{
				alpha *= 1.5 ;
				int k = 0 ;
				likelihood = 0.0 ; 
				for( int i = 0 ; i < _Space ; i++ )
				{
					if( object[i] < 0.0 && ( alpha * object[i] + ws.buf[i] ) < 0.0 ) k++ ;
					temp = ws.eimg[i] + alpha * rat[i] ;
					if ( temp > EMSepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
				}
				if( _CheckStatus ) _EMprintAcceleration( alpha, likelihood, k ) ;       		
				if( k > 0 ) 
				{
					alpha /= 1.5 ;
					continue_loop_1 = false ;
				}
			}
			continue_loop_2 = false ;
		}
		if( continue_loop_2 )
		{
			alpha = alpha_new ;
			temp1 = 0.0 ;
			temp2 = 0.0 ;
			for( int i = 0 ; i < _Space ; i++ ) 
			{
				temp = ws.eimg[i] + alpha * rat[i] ; 
				if( fabs( temp ) > EMSepsilon )
				{ 
					temp1 += ( image[i] * rat[i] / temp - rat[i] ) ;
					temp2 += ( image[i] * rat[i] * rat[i] / ( temp * temp ) ) ;
				}
			}
			alpha_new = alpha + temp1 / temp2 ;
		}
	}
	while( fabs( alpha_new - alpha ) > alpha * 0.1 && continue_loop_2 ) ;
	
	if( _CheckStatus ) _EMprintAcceleration( alpha ) ;
             
	for( int i = 0 ; i < _Space ; i++ ) 
	{
		object[i] = ws.buf[i] + alpha * object[i] ;
		if ( object[i] < 0.0 ) object[i] = 0.0 ;
	}
        
	if( _TrackLikelihood )
	{
		likelihood = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			temp = ws.eimg[i] + alpha * rat[i] ;
			if ( temp > EMSepsilon ) likelihood += ( image[i] * log(temp) - temp ) ;
		}
		_Likelihood.push_back( likelihood ) ;
	}
}
