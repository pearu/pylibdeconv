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
 * Filename:  CGdeconvolver.cc
 */


#include <math.h>
#include "CGdeconvolver.h"


/* public functions */

void CGdeconvolver::init( bool IsApplyIR, bool IsApplyNorm, bool IsTrackLike, bool IsTrackMax, bool IsCheckStatus )
{
	_setControlFlags( IsCheckStatus, IsApplyNorm, IsTrackMax, IsTrackLike ) ;
	
	_ApplyIR = IsApplyIR ;
	_CGIRpenalty = -1.0 ;
	
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



void CGdeconvolver::exportCG( const char * filename )
{
	FILE * fp = fopen( filename, "a+" ) ;
	if( fp )
	{
		_exportLWCG( fp ) ;
		
		fprintf( fp, "%d -> Apply intensity regularization iteratively in the deconvolution loop.\n", (int)_ApplyIR ) ; 
		fprintf( fp, "\n" ) ;
			
		if( _CGIRpenalty > 0.0 ) 
		{
			fprintf( fp, "%e -> Penalty value applied for intensity regularization in the deconvolution loop.\n", _CGIRpenalty ) ;
			fprintf( fp, "\n" ) ;
		}
			
		fclose( fp ) ;
	}
	else
	{
		throw ErrnoError( std::string(filename) ) ;
	}
}



void CGdeconvolver::run( int DimX, int DimY, int DimZ, double * cgr, double * cgp, double * object, 
                         CGdws & ws, unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	double cri = 1.0E+37, max_intensity = 0.0, gamma = -1.0, alpha = 0.0, beta = 0.0, temp1, temp2, temp3 ;

	/* initialize running */
	_CGstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _CGprintStatus( 1 ) ;
	_initPSF( ws.size, cgp, ws.psf_re, ws.psf_im, FrequencySupport, ws.otf ) ;
	_initIMG( max_intensity, cgr, object, SpacialSupport ) ;	
	_FFTplanff->execute( cgr, ws.image_re, ws.image_im ) ;	
	if( _CheckStatus ) _CGprintStatus( 2 ) ;
	
	/* start regularization */
	if( _ApplyIR )
	{
		for( int i = 0 ; i < ws.size ; i++ ) 
		{
			ws.cg_re[i] = ws.image_re[i] * ws.image_re[i] + ws.image_im[i] * ws.image_im[i] ;
		}
		_CGIRpenalty = _CGrunRegularization( ws.size, ws.cg_re, ws.otf ) ;
		if( _CheckStatus ) _CGprintStatus( 3 ) ;
	}
	else
	{
		_CGIRpenalty = 0.0 ;
	}			

	/* start conditioning */
	if( _ConditioningIteration > 0 )
	{
		if( _CheckStatus ) _CGprintStatus( 4 ) ;
		for( int i = 0 ; i < _Space ; i++ ) ws.object0[i] = object[i] ;
		_ConditioningValue = _runConditioning( ws.size, ws.object0, object, cgr, cgp, 
		                     ws.image_re, ws.image_im, ws.psf_re, ws.psf_im, ws.otf, SpacialSupport ) ;
		for( int i = 0 ; i < _Space ; i++ ) object[i] = ws.object0[i] ;
		delete [] ws.object0 ;
		ws.object0 = NULL ;
		if( _CheckStatus ) _CGprintStatus( 5 ) ;
	}
	
	/* initialize arrays in deconvolution loop */
	if( !_TrackLikelihood )
	{
		for( int i = 0 ; i < ws.size ; i++ )
		{
			         temp1 = ws.otf[i] + _ConditioningValue ;
			         temp2 = sqrt( temp1 ) ;
			         temp3 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_im[i] = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
		   	ws.image_re[i] = temp3 ;
		   	  ws.psf_re[i] = ws.psf_re[i] / temp2 ;
		   	  ws.psf_im[i] = ws.psf_im[i] / temp2 ;
		   	     ws.otf[i] = ws.otf[i] / temp1 ;
		}
	}

	/* deconvolution loop */
	if( _CheckStatus ) _CGprintStatus( 6 ) ;	
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _CGprintStatus( 7 ) ;			
		_FFTplanff->execute( object, ws.cg_re, ws.cg_im  ) ;		
		if( _TrackLikelihood )
		{
			_Likelihood.push_back( _getLikelihood( ws.size, ws.image_re, ws.image_im, 
			                       ws.psf_re, ws.psf_im, ws.cg_re, ws.cg_im ) ) ;	
			_CGupdate1( gamma, alpha, beta, cgr, cgp, ws ) ; 
		}
		else
		{
			_CGupdate2( gamma, alpha, beta, cgr, cgp, ws ) ;
		}     		
		for( int i = 0 ; i < _Space ; i++ )
		{
			       cgr[i] = object[i] ;
			    object[i] = object[i] + alpha * cgp[i] ;
			if( object[i] < 0.0 )
			{
				 object[i] = 0.0 ;
				ws.sign[i] = (unsigned char) 0 ;
			} 
			else    ws.sign[i] = (unsigned char) 1 ;
		}
		_getUpdate( object, cgr, SpacialSupport ) ;	
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}
		if( _CheckStatus ) _CGprintStatus( 8 ) ;
	}	
	
	/* end deconvolution */
	_CGfinishRun( ws ) ;	
}



void CGdeconvolver::run( int DimX, int DimY, int DimZ, float * cgr, float * cgp, float * object, 
                         CGsws & ws, unsigned char * SpacialSupport, unsigned char * FrequencySupport )
{
	float  max_intensity = 0.0, gamma = -1.0, alpha = 0.0, beta = 0.0, temp1, temp2, temp3 ;
	double cri = 1.0E+37 ;

	/* initialize running */
	_CGstartRun( DimX, DimY, DimZ, ws ) ;

	/* start initialization */
	if( _CheckStatus ) _CGprintStatus( 1 ) ;
	_initPSF( ws.size, cgp, ws.psf_re, ws.psf_im, FrequencySupport, ws.otf ) ;
	_initIMG( max_intensity, cgr, object, SpacialSupport ) ;
	_FFTplanff->execute( cgr, ws.image_re, ws.image_im ) ;	
	if( _CheckStatus ) _CGprintStatus( 2 ) ;	
	
	/* start regularization */
	if( _ApplyIR )
	{
		for( int i = 0 ; i < ws.size ; i++ )
		{
			ws.cg_re[i] = ws.image_re[i] * ws.image_re[i] + ws.image_im[i] * ws.image_im[i] ;
		}
		_CGIRpenalty = _CGrunRegularization( ws.size, ws.cg_re, ws.otf ) ;
		if( _CheckStatus ) _CGprintStatus( 3 ) ;
	}
	else
	{
		_CGIRpenalty = 0.0 ;
	}			

	/* start conditioning */
	if( _ConditioningIteration > 0 )
	{
		if( _CheckStatus ) _CGprintStatus( 4 ) ;
		for( int i = 0 ; i < _Space ; i++ ) ws.object0[i] = object[i] ;
		_ConditioningValue = _runConditioning( ws.size, ws.object0, object, cgr, cgp, 
		                     ws.image_re, ws.image_im, ws.psf_re, ws.psf_im, ws.otf, SpacialSupport ) ;
		for( int i = 0 ; i < _Space ; i++ ) object[i] = ws.object0[i] ;
		delete [] ws.object0 ;
		ws.object0 = NULL ;
		if( _CheckStatus ) _CGprintStatus( 5 ) ;
	}
	
	/* initialize arrays in deconvolution loop */
	if( !_TrackLikelihood )
	{
		for( int i = 0 ; i < ws.size ; i++ )
		{
			         temp1 = ws.otf[i] + _ConditioningValue ;
			         temp2 = sqrt( temp1 ) ;
			         temp3 = ( ws.image_re[i]*ws.psf_re[i] + ws.image_im[i]*ws.psf_im[i] ) / temp1 ;
			ws.image_im[i] = ( ws.image_im[i]*ws.psf_re[i] - ws.image_re[i]*ws.psf_im[i] ) / temp1 ;
		   	ws.image_re[i] = temp3 ;
		   	  ws.psf_re[i] = ws.psf_re[i] / temp2 ;
		   	  ws.psf_im[i] = ws.psf_im[i] / temp2 ;
		   	     ws.otf[i] = ws.otf[i] / temp1 ;
		}
	}	

	/* deconvolution loop */
	if( _CheckStatus ) _CGprintStatus( 6 ) ;	
	while( _Update.size() < _MaxRunIteration && cri > _Criterion )
	{
		if( _CheckStatus ) _CGprintStatus( 7 ) ;			
		_FFTplanff->execute( object, ws.cg_re, ws.cg_im  ) ;		
		if( _TrackLikelihood )
		{
			_Likelihood.push_back( _getLikelihood( ws.size, ws.image_re, ws.image_im, 
			                       ws.psf_re, ws.psf_im, ws.cg_re, ws.cg_im ) ) ;	
			_CGupdate1( gamma, alpha, beta, cgr, cgp, ws ) ; 
		}
		else
		{
			_CGupdate2( gamma, alpha, beta, cgr, cgp, ws ) ;
		}     		
		for( int i = 0 ; i < _Space ; i++ )
		{
			       cgr[i] = object[i] ;
			    object[i] = object[i] + alpha * cgp[i] ;
			if( object[i] < 0.0 )
			{
				 object[i] = 0.0 ;
				ws.sign[i] = (unsigned char) 0 ;
			} 
			else    ws.sign[i] = (unsigned char) 1 ;
		}
		_getUpdate( object, cgr, SpacialSupport ) ;     		
		if( _Update.size() >= 10 )
		{
			cri = 0.0 ;
			for( int i = 1 ; i <= 10 ; i++ ) cri += _Update[ _Update.size()-i ] ;
			cri /= 10.0 ;
		}		   		   		
		if( _CheckStatus ) _CGprintStatus( 8 ) ;
	}	
	
	/* end deconvolution */
	_CGfinishRun( ws ) ;	
}



/* private functions */

void CGdeconvolver::_CGprintStatus( int stage )
{
        switch( stage )
        {	
		case 1:
			time( &_t0 ) ;
			std::cout << " CGdeconvolver::run starts initialization ... \n" ;
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
			std::cout << " CGdeconvolver::run completes initialization, elapsed "
			          << difftime( _t1, _t0 ) << " seconds.\n" ;
			if( _ApplySpacialSupport )
			{
				std::cout << " CGdeconvolver::run applys spacial support on the deconvolved object iteratively.\n" ;
			}
			break ;

		case 3:
			std::cout << " CGdeconvolver::run applys intensity regularization with penalty = " 
			          << _CGIRpenalty << " iteratively.\n" ;
			break ;

		case 4:
			time( &_t0 ) ;
			std::cout << " CGdeconvolver::run starts conditioning ... \n" ;
			break ;

		case 5:
			time( &_t1 ) ;
			std::cout << " CGdeconvolver::run completes conditioning, elapsed "
			          << difftime( _t1, _t0 ) << " seconds.\n" ;
			break ;

		case 6:
			std::cout << " CGdeconvolver::run starts the loop with conditioning value = " 
			          << _ConditioningValue << ".\n" ;
			break ;

		case 7:
			time( &_t0 ) ;
			break ;
		
                case 8:
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



void CGdeconvolver::_CGstartRun( int DimX, int DimY, int DimZ, CGdws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double)_Space * 3.0 ;
		
	time( &_StartRunTime ) ;
	std::cout << " CGdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " CGdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " CGdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " CGdeconvolution stopping criterion : " << _Criterion << "\n" ;
        	
	time( &_t0 ) ;
	std::cout << " CGdeconvolver::run starts creating FFT plans ... \n" ;
        	
	FFTW3_FFT::Ptr tempff( new FFTW3_FFT( _DimX, _DimY, _DimZ, true,  true, 3 ) ) ;
	_FFTplanff = tempff ;
	
	FFTW3_FFT::Ptr tempbb( new FFTW3_FFT( _DimX, _DimY, _DimZ, false, true, 3 ) ) ;
	_FFTplanbb = tempbb ;
        	
	if( _ConditioningIteration > 0 )
	{
		FFTW3_FFT::Ptr tempf( new FFTW3_FFT( _DimX, _DimY, _DimZ, true,  true, 1 ) ) ;
		_FFTplanf = tempf ;
	
		FFTW3_FFT::Ptr tempb( new FFTW3_FFT( _DimX, _DimY, _DimZ, false, true, 2 ) ) ;
		_FFTplanb = tempb ;
	}
		
	time( &_t1 ) ;
	std::cout << " CGdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;
	
	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 7.0 + ((double)_Space) / 8.0 ) ;
	
	ws.psf_re   = new double[ ws.size ] ;
	ws.psf_im   = new double[ ws.size ] ;
	ws.image_re = new double[ ws.size ] ;
	ws.image_im = new double[ ws.size ] ;
	ws.otf      = new double[ ws.size ] ;
	ws.cg_re    = new double[ ws.size ] ;
	ws.cg_im    = new double[ ws.size ] ;
	ws.sign     = new unsigned char[ _Space ] ;
	if( _ConditioningIteration > 0 )
	{
		ws.object0 = new double[ _Space ] ;
		memory += ( (double)_Space ) ;
	}
	else    ws.object0 = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;
	
	std::cout << " CGdeconvolver::run are using " << (memory/1024.0/128.0) << " Mbytes memory.\n" ;
}



void CGdeconvolver::_CGstartRun( int DimX, int DimY, int DimZ, CGsws & ws )
{
	_setDimensions( DimX, DimY, DimZ ) ;
	double memory = (double)_Space * 3 ;
		
	time( &_StartRunTime ) ;
	std::cout << " CGdeconvolution starts running at " << ctime( &_StartRunTime ) ;
	std::cout << " CGdeconvolution size : " << _DimX << " x " << _DimY << " x " << _DimZ << "\n" ;
	std::cout << " CGdeconvolution max allowed iterations : " << _MaxRunIteration << "\n" ;
	std::cout << " CGdeconvolution stopping criterion : " << _Criterion << "\n" ;
        	
	time( &_t0 ) ;
	std::cout << " CGdeconvolver::run starts creating FFT plans ... \n" ;
        	
	FFTW3_FFT::Ptr tempff( new FFTW3_FFT( _DimX, _DimY, _DimZ, true,  false, 3 ) ) ;
	_FFTplanff = tempff ;
	
	FFTW3_FFT::Ptr tempbb( new FFTW3_FFT( _DimX, _DimY, _DimZ, false, false, 3 ) ) ;
	_FFTplanbb = tempbb ;
		
	if( _ConditioningIteration > 0 )
	{
		FFTW3_FFT::Ptr tempf( new FFTW3_FFT( _DimX, _DimY, _DimZ, true,  false, 1 ) ) ;
		_FFTplanf = tempf ;
	
		FFTW3_FFT::Ptr tempb( new FFTW3_FFT( _DimX, _DimY, _DimZ, false, false, 2 ) ) ;
		_FFTplanb = tempb ;
	}
		
	time( &_t1 ) ;
	std::cout << " CGdeconvolver::run completes creating FFT plans, elapsed "
	          << difftime( _t1, _t0 ) << " seconds.\n" ;
	
	ws.size = _FFTplanf->FFTsize() ;
	memory += ( (double)ws.size * 7.0 + ((double)_Space) / 8.0 ) ;
	
	ws.psf_re   = new float[ ws.size ] ;
	ws.psf_im   = new float[ ws.size ] ;
	ws.image_re = new float[ ws.size ] ;
	ws.image_im = new float[ ws.size ] ;
	ws.otf      = new float[ ws.size ] ;
	ws.cg_re    = new float[ ws.size ] ;
	ws.cg_im    = new float[ ws.size ] ;
	ws.sign     = new unsigned char[ _Space ] ;
	if( _ConditioningIteration > 0 )
	{
		ws.object0 = new float[ _Space ] ;
		memory += ( (double)_Space ) ;
	}
	else    ws.object0 = NULL ;
	
	if( _Update.size()  > 0 ) _Update.clear() ;
	if( _Likelihood.size() > 0 ) _Likelihood.clear() ;
	if( _ObjectMax.size()  > 0 ) _ObjectMax.clear() ;
	
	std::cout << " CGdeconvolver::run are using " << (memory/1024.0/256.0) << " Mbytes memory.\n" ;
}



void CGdeconvolver::_CGfinishRun( CGdws & ws )
{
	if( ws.psf_re   != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im   != NULL ) delete [] ws.psf_im ;
	if( ws.image_re != NULL ) delete [] ws.image_re ;
	if( ws.image_im != NULL ) delete [] ws.image_im ;
	if( ws.otf      != NULL ) delete [] ws.otf ;
	if( ws.cg_re    != NULL ) delete [] ws.cg_re ;
	if( ws.cg_im    != NULL ) delete [] ws.cg_im ;
	if( ws.sign     != NULL ) delete [] ws.sign ;
		
	time( &_StopRunTime ) ;
	std::cout << " CGdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}  



void CGdeconvolver::_CGfinishRun( CGsws & ws )
{
	if( ws.psf_re   != NULL ) delete [] ws.psf_re ;
	if( ws.psf_im   != NULL ) delete [] ws.psf_im ;
	if( ws.image_re != NULL ) delete [] ws.image_re ;
	if( ws.image_im != NULL ) delete [] ws.image_im ;
	if( ws.otf      != NULL ) delete [] ws.otf ;
	if( ws.cg_re    != NULL ) delete [] ws.cg_re ;
	if( ws.cg_im    != NULL ) delete [] ws.cg_im ;
	if( ws.sign     != NULL ) delete [] ws.sign ;
		
	time( &_StopRunTime ) ;
	std::cout << " CGdeconvolution finish running at " << ctime( &_StopRunTime ) ;
}  



void CGdeconvolver::_CGupdate1( double & gamma, double & alpha, double & beta, double * cgr, double * cgp, CGdws & ws )
{
	double temp1, temp2, temp3, temp4 ;
		
	for( int i = 0 ; i < ws.size ; i++ )
	{
		      temp1 = ws.otf[i] + _ConditioningValue ;
		      temp2 = ( ws.image_re[i] * ws.psf_re[i] + ws.image_im[i] * ws.psf_im[i] ) / temp1 ;
		      temp3 = ( ws.image_im[i] * ws.psf_re[i] - ws.image_re[i] * ws.psf_im[i] ) / temp1 ;
		      temp4 = ws.otf[i] / temp1 ;
		ws.cg_re[i] = temp2 - ( temp4 + _CGIRpenalty ) * ws.cg_re[i] ;
		ws.cg_im[i] = temp3 - ( temp4 + _CGIRpenalty ) * ws.cg_im[i] ; 
	}
	
 	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
 	
 	if( gamma < 0.0 ) 
 	{
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) 
 		{
 			    cgp[i] = cgr[i] ;
 			ws.sign[i] = ( unsigned char ) 1 ;
 			   gamma  += cgr[i] * cgr[i] ;
 		}
 	}
 	else
 	{
 		temp1 = gamma ;
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) gamma += cgr[i] * cgr[i] ;
 		beta = gamma / temp1 ;
 		for( int i = 0 ; i < _Space ; i++ ) cgp[i] = cgr[i] + beta * cgp[i] ;
 		_FFTplanff->execute( cgp, ws.cg_re, ws.cg_im ) ;
 	}
 	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = sqrt( ws.otf[i] + _ConditioningValue ) ; 
		       temp2 = ws.psf_re[i] / temp1 ;
		       temp3 = ws.psf_im[i] / temp1 ;
		       temp4 = ws.cg_re[i] * temp3 + ws.cg_im[i] * temp2 ;
		 ws.cg_re[i] = ws.cg_re[i] * temp2 - ws.cg_im[i] * temp3 ;
		 ws.cg_im[i] = temp4 ;
	}
		
	temp1 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp1 += ( cgr[i] * cgp[i] * ((double) ws.sign[i]) ) ;
     	
	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
	temp2 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp2 += ( ((double) ws.sign[i]) * cgr[i] * cgr[i] ) ;
	if( _ApplyIR )
	{
		temp3 = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ ) temp3 += ( ((double) ws.sign[i]) * cgp[i] * cgp[i] ) ;
		temp2 += temp3 * _CGIRpenalty ;
	}
     	
	alpha = temp1 / temp2 ;
}



void CGdeconvolver::_CGupdate1( float & gamma, float & alpha, float & beta, float * cgr, float * cgp, CGsws & ws )
{
	float temp1, temp2, temp3, temp4 ;
		
	for( int i = 0 ; i < ws.size ; i++ )
	{
		      temp1 = ws.otf[i] + _ConditioningValue ;
		      temp2 = ( ws.image_re[i] * ws.psf_re[i] + ws.image_im[i] * ws.psf_im[i] ) / temp1 ;
		      temp3 = ( ws.image_im[i] * ws.psf_re[i] - ws.image_re[i] * ws.psf_im[i] ) / temp1 ;
		      temp4 = ws.otf[i] / temp1 ;
		ws.cg_re[i] = temp2 - ( temp4 + _CGIRpenalty ) * ws.cg_re[i] ;
		ws.cg_im[i] = temp3 - ( temp4 + _CGIRpenalty ) * ws.cg_im[i] ; 
	}
	
 	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
 	
 	if( gamma < 0.0 ) 
 	{
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) 
 		{
 			    cgp[i] = cgr[i] ;
 			ws.sign[i] = ( unsigned char ) 1 ;
 			   gamma  += cgr[i] * cgr[i] ;
 		}
 	}
 	else
 	{
 		temp1 = gamma ;
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) gamma += cgr[i] * cgr[i] ;
 		beta = gamma / temp1 ;
 		for( int i = 0 ; i < _Space ; i++ ) cgp[i] = cgr[i] + beta * cgp[i] ;
 		_FFTplanff->execute( cgp, ws.cg_re, ws.cg_im ) ;
 	}
 	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = sqrt( ws.otf[i] + _ConditioningValue ) ; 
		       temp2 = ws.psf_re[i] / temp1 ;
		       temp3 = ws.psf_im[i] / temp1 ;
		       temp4 = ws.cg_re[i] * temp3 + ws.cg_im[i] * temp2 ;
		 ws.cg_re[i] = ws.cg_re[i] * temp2 - ws.cg_im[i] * temp3 ;
		 ws.cg_im[i] = temp4 ;
	}
		
	temp1 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp1 += ( cgr[i] * cgp[i] * ((float) ws.sign[i]) ) ;
     	
	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
	temp2 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp2 += ( ((float) ws.sign[i]) * cgr[i] * cgr[i] ) ;
	if( _ApplyIR )
	{
		temp3 = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ ) temp3 += ( ((float) ws.sign[i]) * cgp[i] * cgp[i] ) ;
		temp2 += temp3 * _CGIRpenalty ;
	}
     	
	alpha = temp1 / temp2 ;
}



void CGdeconvolver::_CGupdate2( double & gamma, double & alpha, double & beta, double * cgr, double * cgp, CGdws & ws )
{
	double temp1, temp2, temp3 ;
		
	for( int i = 0 ; i < ws.size ; i++ )
	{
		ws.cg_re[i] = ws.image_re[i] - ( ws.otf[i] + _CGIRpenalty ) * ws.cg_re[i] ;
		ws.cg_im[i] = ws.image_im[i] - ( ws.otf[i] + _CGIRpenalty ) * ws.cg_im[i] ; 
	}
	
 	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
 	
 	if( gamma < 0.0 ) 
 	{
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) 
 		{
 			    cgp[i] = cgr[i] ;
 			ws.sign[i] = ( unsigned char ) 1 ;
 			   gamma  += ( cgr[i] * cgr[i] ) ;
 		}
 	}
 	else
 	{
 		temp1 = gamma ;
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) gamma += cgr[i] * cgr[i] ;
 		beta = gamma / temp1 ;
 		for( int i = 0 ; i < _Space ; i++ ) cgp[i] = cgr[i] + beta * cgp[i] ;
 		_FFTplanff->execute( cgp, ws.cg_re, ws.cg_im ) ;
 	}
 	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = ws.cg_re[i] * ws.psf_re[i] - ws.cg_im[i] * ws.psf_im[i] ;
		 ws.cg_im[i] = ws.cg_re[i] * ws.psf_im[i] + ws.cg_im[i] * ws.psf_re[i] ;
		 ws.cg_re[i] = temp1 ;
	}
		
	temp1 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp1 += ( cgr[i] * cgp[i] * ((double) ws.sign[i]) ) ;
     	
	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
	temp2 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp2 += ( ((double) ws.sign[i]) * cgr[i] * cgr[i] ) ; 
	if( _ApplyIR )
	{
		temp3 = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ ) temp3 += ( ((float) ws.sign[i]) * cgp[i] * cgp[i] ) ;
		temp2 += temp3 * _CGIRpenalty ;	
	}
     		
	alpha = temp1 / temp2 ;
}



void CGdeconvolver::_CGupdate2( float & gamma, float & alpha, float & beta, float * cgr, float * cgp, CGsws & ws )
{
	float temp1, temp2, temp3 ;
		
	for( int i = 0 ; i < ws.size ; i++ )
	{
		ws.cg_re[i] = ws.image_re[i] - ( ws.otf[i] + _CGIRpenalty ) * ws.cg_re[i] ;
		ws.cg_im[i] = ws.image_im[i] - ( ws.otf[i] + _CGIRpenalty ) * ws.cg_im[i] ; 
	}
	
 	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
 	
 	if( gamma < 0.0 ) 
 	{
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) 
 		{
 			    cgp[i] = cgr[i] ;
 			ws.sign[i] = ( unsigned char ) 1 ;
 			   gamma  += ( cgr[i] * cgr[i] ) ;
 		}
 	}
 	else
 	{
 		temp1 = gamma ;
 		gamma = 0.0 ;
 		for( int i = 0 ; i < _Space ; i++ ) gamma += cgr[i] * cgr[i] ;
 		beta = gamma / temp1 ;
 		for( int i = 0 ; i < _Space ; i++ ) cgp[i] = cgr[i] + beta * cgp[i] ;
 		_FFTplanff->execute( cgp, ws.cg_re, ws.cg_im ) ;
 	}
 	
	for( int i = 0 ; i < ws.size ; i++ )
	{
		       temp1 = ws.cg_re[i] * ws.psf_re[i] - ws.cg_im[i] * ws.psf_im[i] ;
		 ws.cg_im[i] = ws.cg_re[i] * ws.psf_im[i] + ws.cg_im[i] * ws.psf_re[i] ;
		 ws.cg_re[i] = temp1 ;
	}
		
	temp1 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp1 += ( cgr[i] * cgp[i] * ((float) ws.sign[i]) ) ;
     	
	_FFTplanbb->execute( ws.cg_re, ws.cg_im, cgr ) ;
	temp2 = 0.0 ;
	for( int i = 0 ; i < _Space ; i++ ) temp2 += ( ((float) ws.sign[i]) * cgr[i] * cgr[i] ) ; 
	if( _ApplyIR )
	{
		temp3 = 0.0 ;
		for( int i = 0 ; i < _Space ; i++ ) temp3 += ( ((float) ws.sign[i]) * cgp[i] * cgp[i] ) ;
		temp2 += temp3 * _CGIRpenalty ;	
	}
     		
	alpha = temp1 / temp2 ;
}



double CGdeconvolver::_CGrunRegularization( int size, double * img, double * otf )
{
	double R   = 0.61803399 ;
	double C   = 1.0 - R ;
	double tol = 1.0E-10 ;
	double gcv = 1.0E+37 ;
	double last_gcv = gcv ;
	double x0, x1, x2, x3, f1, f2, gcv1, gcv2 ;

	x0 = 1.0 ;
	while( gcv <= last_gcv )
	{
		last_gcv = gcv ;
		x0 = x0 * 0.1 ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x0 * x0 * img[i] / ( otf[i] + x0 ) / ( otf[i] + x0 ) ) ;
			gcv2 += ( x0 / ( otf[i] + x0 ) ) ;
		}
		gcv = gcv1 / gcv2 / gcv2 ;
	}        	       	       

	x1 = x0 * 10.0 ;
	x3 = x0 * 100.0 ;    	
	if( fabs(x3-x1) > fabs(x1-x0) )
	{
		f1 = last_gcv ;
		x2 = x1 + C * ( x3 - x1 ) ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x2 * x2 * img[i] / ( otf[i] + x2 ) / ( otf[i] + x2 ) ) ;
			gcv2 += ( x2 / ( otf[i] + x2 ) ) ;
		}
		f2 = gcv1 / gcv2 / gcv2 ;
	}
	else
	{
		x2 = x1 ;
		f2 = last_gcv ;
		x1 = x2 - C * ( x2 - x0 ) ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x1 * x1 * img[i] / ( otf[i] + x1 ) / ( otf[i] + x1 ) ) ;
			gcv2 += ( x1 / ( otf[i] + x1 ) ) ;
		}
		f1 = gcv1 / gcv2 / gcv2 ;
	}
       	
	while( fabs(x3-x0) > tol * ( fabs(x1) + fabs(x2) ) )
	{
		if ( f2 < f1 )
		{
			x0 = x1 ;
			x1 = x2 ;
			x2 = R * x1 + C * x3 ;
			f1 = f2 ;
			gcv1 = 0.0 ;
			gcv2 = 0.0 ;
			for( int i = 0 ; i < size ; i++ )
			{ 
				gcv1 += ( x2 * x2 * img[i] / ( otf[i] + x2 ) / ( otf[i] + x2 ) ) ;
				gcv2 += ( x2 / ( otf[i] + x2 ) ) ;
			}
			f2 = gcv1 / gcv2 / gcv2 ;
		}
		else
		{ 
			x3 = x2 ; 
			x2 = x1 ;
			x1 = R * x2 + C * x0 ;
			f2 = f1 ;
			gcv1 = 0.0 ;
			gcv2 = 0.0 ;
			for( int i = 0 ; i < size ; i++ )
			{ 
				gcv1 += ( x1 * x1 * img[i] / ( otf[i] + x1 ) / ( otf[i] + x1 ) ) ;
				gcv2 += ( x1 / ( otf[i] + x1 ) ) ;
			}
			f1 = gcv1 / gcv2 / gcv2 ;
		}
	}

	if( f1 < f2 )
	{
		return x1 ;
	}
	else 
	{
		return x2 ;
	}
}



double CGdeconvolver::_CGrunRegularization( int size, float * img, float * otf )
{
	double R   = 0.61803399 ;
	double C   = 1.0 - R ;
	double tol = 1.0E-10 ;
	double gcv = 1.0E+37 ;
	double last_gcv = gcv ;
	double x0, x1, x2, x3, f1, f2, gcv1, gcv2 ;

	x0 = 1.0 ;
	while( gcv <= last_gcv )
	{
		last_gcv = gcv ;
		x0 = x0 * 0.1 ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x0 * x0 * img[i] / ( otf[i] + x0 ) / ( otf[i] + x0 ) ) ;
			gcv2 += ( x0 / ( otf[i] + x0 ) ) ;
		}
		gcv = gcv1 / gcv2 / gcv2 ;
	}        	       	       

	x1 = x0 * 10.0 ;
	x3 = x0 * 100.0 ;    	
	if( fabs(x3-x1) > fabs(x1-x0) )
	{
		f1 = last_gcv ;
		x2 = x1 + C * ( x3 - x1 ) ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x2 * x2 * img[i] / ( otf[i] + x2 ) / ( otf[i] + x2 ) ) ;
			gcv2 += ( x2 / ( otf[i] + x2 ) ) ;
		}
		f2 = gcv1 / gcv2 / gcv2 ;
	}
	else
	{
		x2 = x1 ;
		f2 = last_gcv ;
		x1 = x2 - C * ( x2 - x0 ) ;
		gcv1 = 0.0 ;
		gcv2 = 0.0 ;
		for( int i = 0 ; i < size ; i++ )
		{ 
			gcv1 += ( x1 * x1 * img[i] / ( otf[i] + x1 ) / ( otf[i] + x1 ) ) ;
			gcv2 += ( x1 / ( otf[i] + x1 ) ) ;
		}
		f1 = gcv1 / gcv2 / gcv2 ;
	}
       	
	while( fabs(x3-x0) > tol * ( fabs(x1) + fabs(x2) ) )
	{
		if ( f2 < f1 )
		{
			x0 = x1 ;
			x1 = x2 ;
			x2 = R * x1 + C * x3 ;
			f1 = f2 ;
			gcv1 = 0.0 ;
			gcv2 = 0.0 ;
			for( int i = 0 ; i < size ; i++ )
			{ 
				gcv1 += ( x2 * x2 * img[i] / ( otf[i] + x2 ) / ( otf[i] + x2 ) ) ;
				gcv2 += ( x2 / ( otf[i] + x2 ) ) ;
			}
			f2 = gcv1 / gcv2 / gcv2 ;
		}
		else
		{ 
			x3 = x2 ; 
			x2 = x1 ;
			x1 = R * x2 + C * x0 ;
			f2 = f1 ;
			gcv1 = 0.0 ;
			gcv2 = 0.0 ;
			for( int i = 0 ; i < size ; i++ )
			{ 
				gcv1 += ( x1 * x1 * img[i] / ( otf[i] + x1 ) / ( otf[i] + x1 ) ) ;
				gcv2 += ( x1 / ( otf[i] + x1 ) ) ;
			}
			f1 = gcv1 / gcv2 / gcv2 ;
		}
	}

	if( f1 < f2 )
	{
		return x1 ;
	}
	else 
	{
		return x2 ;
	}
}
