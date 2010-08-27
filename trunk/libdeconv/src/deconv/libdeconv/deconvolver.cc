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
 * Filename:  deconvolver.cc
 */
 

#include "deconvolver.h"
#include "FFTW3fft.h"
 
 
/* public functions */

unsigned int deconvolver::exportUpdateTrack( double * vtrack ) 
{
	if( _Update.size() > 0 )
	{
		for( unsigned int i = 0 ; i < _Update.size() ; i++ ) vtrack[i] = _Update[i] ;
	}
	
	return ( _Update.size() ) ;
}



unsigned int deconvolver::exportObjectMaxTrack( double * vtrack )
{
	if( _ObjectMax.size() > 0 ) 
	{
		for( unsigned int i = 0 ; i < _ObjectMax.size() ; i++ ) vtrack[i] = _ObjectMax[i] ;
	}
	
	return ( _ObjectMax.size() ) ;
}



unsigned int deconvolver::exportLikelihoodTrack( double * vtrack ) 
{
	if( _Likelihood.size() > 0 )
	{
		for( unsigned int i = 0 ; i < _Likelihood.size() ; i++ ) vtrack[i] = _Likelihood[i] ;
	}
	
	return ( _Likelihood.size() ) ;
}



/* protected functions */

bool deconvolver::_IsPowerOf2( int num )
{
	int newnum, pcnt ;
        
	pcnt = 0 ;
	newnum = num ;
        
	while( newnum > 0 ) 
	{
		newnum = newnum >> 1 ; 
		pcnt++ ; 
	}
	
	pcnt-- ;
	newnum = 1 << pcnt;

	if( newnum == num ) 
	{
		return true ;
	}
	else    return false ;
}



void deconvolver::_exportCommon( FILE * fp ) 
{
	if( _StartRunTime > 0 && _StopRunTime > 0 ) 
	{
		fprintf( fp, "deconvolver::run() began at %s", ctime( &_StartRunTime ) ) ;
		fprintf( fp, "deconvolver::run() ended at %s", ctime( &_StopRunTime ) ) ;
		fprintf( fp, "\n" ) ;		
		fprintf( fp, "%d -> Frequency_Support was applied on the FFT of the input PSF.\n", ((int) _ApplyFrequencySupport) ) ; 
		fprintf( fp, "%d -> Spacial_Support was applied on the devonvolved object iteratively.\n", ((int) _ApplySpacialSupport) ) ;
		fprintf( fp, "\n" ) ;
	}
		
	fprintf( fp, "%d -> Dimension X of deconvolution size\n", _DimX ) ; 
	fprintf( fp, "%d -> Dimension Y of deconvolution size\n", _DimY ) ;
	fprintf( fp, "%d -> Dimension Z of deconvolution size\n", _DimZ ) ;
	fprintf( fp, "\n" ) ;

	fprintf( fp, "%d -> Max_Allowed_Iterations in the deconvolution loop\n", _MaxRunIteration ) ;
	fprintf( fp, "%d -> Actually_Run_Iterations in the deconvolution loop\n", (int)_Update.size() ) ;
	fprintf( fp, "%e -> Stop_Criterion_to_Terminate the deconvolution loop\n",  _Criterion ) ;
	fprintf( fp, "\n" ) ;
	
	fprintf( fp, "%d -> Apply Normalization on the input image and deconvolved object.\n", ((int) _ApplyNormalization) ) ;
	fprintf( fp, "%d -> Track Deconvolved_Object_Max_Value iteratively in deconvolution loop.\n", ((int) _TrackMaxInObject) ) ;
	fprintf( fp, "%d -> Track Likelihood_Value iteratively in deconvolution loop.\n", ((int) _TrackLikelihood) ) ;
	fprintf( fp, "%d -> Check deconvolution_Process_Running_Status in Terminal.\n", ((int) _CheckStatus) ) ;
	fprintf( fp, "\n" ) ;
	
}

 

void deconvolver::_setDimensions( int DimX, int DimY, int DimZ )
{
	if( DimX > 0 && _IsPowerOf2(DimX) && DimY > 0 && _IsPowerOf2(DimY) && DimZ > 0 && _IsPowerOf2(DimZ)  )
	{
		_DimX  = DimX ;
		_DimY  = DimY ;
		_DimZ  = DimZ ;
		_Space = DimX * DimY * DimZ ;
	} 
	else	throw DimensionError( DimX, DimY, DimZ ) ;
}



void deconvolver::_setControlFlags(  bool IsCheck, bool IsApply, bool IsTrackMax, bool IsTrackLike )
{
	_CheckStatus        = IsCheck ;
	_TrackLikelihood    = IsTrackLike ;
	_TrackMaxInObject   = IsTrackMax ;	
	_ApplyNormalization = IsApply ;
	if( IsApply ) _TrackMaxInObject = true ;
	_ApplySpacialSupport = false ;
	_ApplyFrequencySupport = false ;
}



void deconvolver::_initPSF( int size, double * psf, double * psf_re, double * psf_im, unsigned char * FrequencySupport, double * otf )
{
	fft3d( _DimX, _DimY, _DimZ, psf, psf_re, psf_im ) ;
	
	if( FrequencySupport != NULL )
	{
		_ApplyFrequencySupport = true ;
		for( int i = 0 ; i < size ; i++ )
		{
			psf_re[i] *= ((double) FrequencySupport[i]) ;
			psf_im[i] *= ((double) FrequencySupport[i]) ;
		}
	}
	
	if( otf != NULL )
	{
		for( int i = 0 ; i < size ; i++ ) 
		{
			otf[i] = psf_re[i] * psf_re[i] + psf_im[i] * psf_im[i] ;
		}
	}
}



void deconvolver::_initPSF( int size, float * psf, float * psf_re, float * psf_im, unsigned char * FrequencySupport, float * otf )
{
	fft3d( _DimX, _DimY, _DimZ, psf, psf_re, psf_im ) ;
	
	if( FrequencySupport != NULL )
	{
		_ApplyFrequencySupport = true ;
		for( int i = 0 ; i < size ; i++ )
		{
			psf_re[i] *= ((float) FrequencySupport[i]) ;
			psf_im[i] *= ((float) FrequencySupport[i]) ;
		}
	}
	
	if( otf != NULL )
	{
		for( int i = 0 ; i < size ; i++ ) 
		{
			otf[i] = psf_re[i] * psf_re[i] + psf_im[i] * psf_im[i] ;
		}
	}	
}



void deconvolver::_initIMG( double & max_intensity, double * image, double * object, unsigned char * SpacialSupport )
{
	if( SpacialSupport != NULL ) _ApplySpacialSupport = true ;
	
	if( _ApplyNormalization )
	{
		max_intensity = object[0] ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] > max_intensity ) max_intensity = object[i] ;
		}
		for( int i = 0 ; i < _Space ; i++ ) object[i] /= max_intensity ;
		
		max_intensity = image[0] ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( image[i] > max_intensity ) max_intensity = image[i] ;
		}
		for( int i = 0 ; i < _Space ; i++ ) image[i] /= max_intensity ;
	}
}



void deconvolver::_initIMG( float & max_intensity, float * image, float * object, unsigned char * SpacialSupport )
{
	if( SpacialSupport != NULL ) _ApplySpacialSupport = true ;
	
	if( _ApplyNormalization )
	{
		max_intensity = object[0] ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] > max_intensity ) max_intensity = object[i] ;
		}
		for( int i = 0 ; i < _Space ; i++ ) object[i] /= max_intensity ;
		
		max_intensity = image[0] ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( image[i] > max_intensity ) max_intensity = image[i] ;
		}
		for( int i = 0 ; i < _Space ; i++ ) image[i] /= max_intensity ;
	}
}



void deconvolver::_getUpdate( double * object, double * last_object, unsigned char * SpacialSupport )
{
	if( SpacialSupport != NULL )
	{
		for( int i = 0 ; i < _Space ; i++ ) object[i] *= ((double) SpacialSupport[i]) ;
	}
     	    		
	if( _TrackMaxInObject )
	{
		double max_intensity = object[0] ;
     		
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] > max_intensity  ) max_intensity = object[i] ;
		}
     		
		_ObjectMax.push_back( max_intensity ) ;

		if( _ApplyNormalization )
		{
			for( int i = 0 ; i < _Space ; i++ ) object[i] /= max_intensity ;
		}
	}
	     			
	double temp1 = 0.0 ;
	double temp2 = 0.0 ;
     	
	for( int i = 0 ; i < _Space ; i++ )
	{
		temp1 += ( object[i] - last_object[i] )* ( object[i] - last_object[i] ) ;
		temp2 += ( object[i] * object[i] ) ;
	}

	_Update.push_back( (temp1/temp2) ) ;
}

       
        
void deconvolver::_getUpdate( float * object, float * last_object, unsigned char * SpacialSupport )
{
	if( SpacialSupport != NULL )
	{
		for( int i = 0 ; i < _Space ; i++ ) object[i] *= ((float) SpacialSupport[i]) ;
	}    		
     	
	if( _TrackMaxInObject )
	{
		float max_intensity = object[0] ;
		for( int i = 0 ; i < _Space ; i++ )
		{
			if( object[i] > max_intensity  ) max_intensity = object[i] ;
		}
		
		_ObjectMax.push_back( max_intensity ) ;
     		
		if( _ApplyNormalization )
		{
			for( int i = 0 ; i < _Space ; i++ ) object[i] /= max_intensity ;
		}
	}     			
     		
	double temp1 = 0.0 ;
	double temp2 = 0.0 ;
     	
	for( int i = 0 ; i < _Space ; i++ )
	{
		temp1 += ( object[i] - last_object[i] )* ( object[i] - last_object[i] ) ;
		temp2 += ( object[i] * object[i] ) ;
	}

	_Update.push_back( (temp1/temp2) ) ;
}
