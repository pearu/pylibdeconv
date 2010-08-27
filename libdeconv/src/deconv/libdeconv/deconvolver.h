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
 * Filename:  deconvolver.h
 */
 
 
#ifndef DECONVOLVER_H
#define DECONVOLVER_H


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include "MYdef.h"
#include "MYerror.h"


/*
 *	==========================================================================================
 *	deconvolver is a base class for the LWdeconvoler, CGdeconvolver and EMdeconvolver classes
 *	since they use some same initial parameters to control how to run a deconvolution process.
 *	==========================================================================================
 *
 *
 *	-------------------------------------
 *	Apply Normailization : <_IsApplyNorm>
 *	-------------------------------------
 *
 *		<_IsApplyNorm>, it is the apply_normalization indicator; its default value is false; 
 *		                If it is true, 
 *		                the input image and iteratively estimated object will be normalized during a deconvolution.
 *
 *	
 *	----------------------------------------
 *	Track Likelihood Value : <_IsTrackLike>
 *	----------------------------------------
 *
 *		<_IsTrackLike>, it is the track_likelihood indicator; its default value is false;
 *		                If it is true, 
 *		                the likelihood value will be tracked iteratively during a deconvolution;
 *		Warnning: be aware that the deconvolution process will be slowed if <_IsTrackLike> is true.
 *
 *
 *	-----------------------------------------------------------------------
 *	Track Max Intensity in an iteratively Estimated Object : <_IsTrackLike>
 *	-----------------------------------------------------------------------
 *
 *		<_IsTrackMax>, it is the track_max_intensity indicator; its default value is false; 
 *		               If it is true, 
 *		               the max intensity in an iteratively estimated object will be tracked.
 *		Warnning: be aware that <_IsTrackMax> will be automatically set to be true if <IsApplyNorm> is true
 *
 *
 *	-------------------------------------------------------------
 *	Track Deconvolution Progress in the Terminal : <_CheckStatus>
 *	-------------------------------------------------------------
 *	
 *		<_CheckStatus>, it is the check_program_running indicator; its default value is true;
 *		                If it is true, 
 *		                the deconvolution progress will be printed in the terminal;
 *
 *
 *	-----------------------------------------------------------------
 *	Criterion to Stop Deconvolution Process  : <_MaxRunIteration>
 *	-----------------------------------------------------------------
 *
 *		<_Criterion>, it is the criterion to stop a deconvolution process; its default value is 1.0e-7;
 *		              Deconvolution process will be stopped when the average value of 
 *		              "||object(k)-object(k-1)||/||object(k)|| over 10 past iterations is less than <_Criterion>.  
 *		              
 *
 *	-----------------------------------------------------------------
 *	Number of Max Allowed Deconvolved Iterations : <_MaxRunIteration>
 *	-----------------------------------------------------------------
 *
 *	<_MaxRunIteration>, it is the number of max allowed deconvolved iterations; its default value is 1000;
 *	                    Deconvolution process will be stopped when <_MaxRunIteration> iteraions have been 
 *	                    executed in the deconvolution main loop even if the criterion is not achieved.
 */

class DimensionError : public Error
{
	public:
	DimensionError( int DimX, int DimY, int DimZ )
	{
		_error << " Deconvolver::run() : each dimension must be positive and power of 2.\n"
		       << " dimensions was set -> DimX = " << DimX 
		       << " , DimY = " << DimY << " , DimZ = " << DimZ << ".\n" ; 
	}
} ;
  
 
		                 
class deconvolver
{
	public:
	virtual ~deconvolver() {}
	
	
	/*
	 *	Get protected members - dimensions of the deconvolution
	 *	DimX()  returns the fastest varying dimension of the cubic image.
	 *	DimY()  returns the middle          dimension of the cubic image.
	 *	DimY()  returns the slowest varying dimension of the cubic image.
	 *	Space() returns the size of the deconvolution and equals to DimX()*DimY()*DimZ().
	 */
	int     DimX()                { return _DimX ;               }
	int     DimY()                { return _DimY ;               }
	int     DimZ()                { return _DimZ ;               }
	int     Space()               { return _Space ;              }

	
	/*	Get protected members - control flags in the deconvolution
	 *	CheckStatus()        returns the check_program_running indicator.
	 *	ApplyNormalization() returns the apply_normalization   indicator.
	 *	TrackMaxInObject()   returns the track_max_intensity   indicator.
	 *	TrackLikelihood()    returns the track_likelihood      indicator.
	 */	
	bool    CheckStatus()         { return _CheckStatus ;        }
	bool    ApplyNormalization()  { return _ApplyNormalization ; }
	bool    TrackMaxInObject()    { return _TrackMaxInObject ;   }
	bool    TrackLikelihood()     { return _TrackLikelihood ;    }
	
	
	/*
	 *	Get protected members - when to stop the deconvolution
	 *	MaxRunIteration() returns the number of max allowed interations for the deconvolution.
	 *	Criterion() 	  returns the criterion to stop the deconvolution.
	 */
	unsigned int  MaxRunIteration()     { return _MaxRunIteration ;    }
	double        Criterion()           { return _Criterion ;          }
	
	
	/*
	 *	Set the number of max allowed deconvolved iterations
	 *	Input:
	 *		iter, it is the number of max allowed iterations and its default value is 1000.
	 */
	void    setMaxRunIteration( unsigned int iter = 1000 ) { _MaxRunIteration = iter ; } 
	
	
	/*
	 *	Set the criterion to stop the deconvolution
	 *	Input:
	 *		cri, it is the criterion and its default value is 1.0e-7.
	 */     
	void    setCriterion( double cri = 1.0e-7 ) { _Criterion = cri ;        }
        
        
	/* 
	 *	Export the update tracking array in deconvolution
	 *	Input:
	 *		vtrack, which points to the array to hold the exported data; 
	 *		        its size must be big enough and the safe size is MaxRunIteration(). 
	 */
	unsigned int  exportUpdateTrack( double * vtrack ) ;
        
        
	/* 
	 *	Export the max intensity tracking array in deconvolution
	 *	Input:
	 *		vtrack, which points to the array to hold the exported data; 
	 *		        its size must be big enough and the safe size is MaxRunIteration(). 
	 */
	unsigned int  exportObjectMaxTrack( double * vtrack ) ;
	
	
	/* 
	 *	Export the likelihood tracking array in deconvolution
	 *	Input:
	 *		vtrack, which points to the array to hold the exported data; 
	 *		        its size must be big enough and the safe size is MaxRunIteration(). 
	 */
	unsigned int  exportLikelihoodTrack( double * vtrack ) ;
	

	protected:
	int                     _DimX ;
	int                     _DimY ;
	int                     _DimZ ;
	int                     _Space ;        
	unsigned int            _MaxRunIteration ;
	double                  _Criterion ;        
	bool                    _CheckStatus ;
	bool                    _ApplyNormalization ;
	bool                    _TrackMaxInObject ;
	bool                    _TrackLikelihood ;
	std::vector< double >   _Update ; 
	std::vector< double >   _ObjectMax ;
	std::vector< double >   _Likelihood ;              
	time_t                  _StartRunTime ;
	time_t                  _StopRunTime ;
	time_t                  _t0, _t1 ;
	bool                    _ApplySpacialSupport ;
	bool                    _ApplyFrequencySupport ;
        
	bool  _IsPowerOf2( int num ) ;
        
	void  _exportCommon( FILE * fp ) ;
        
	void  _setDimensions( int DimX, int DimY, int DimZ ) ;
        
	void  _setControlFlags( bool IsCheck, bool IsApply, bool IsTrackMax, bool IsTrackLike ) ;
	
	void  _initPSF( int size, double * psf, double * psf_re, double * psf_im, unsigned char * FrequencySupport, double * otf = NULL ) ;	 
	void  _initPSF( int size, float  * psf, float  * psf_re, float  * psf_im, unsigned char * FrequencySupport, float  * otf = NULL ) ;
        
	void  _initIMG( double & max_intensity, double * image, double * object, unsigned char * SpacialSupport ) ;
	void  _initIMG( float  & max_intensity, float  * image, float  * object, unsigned char * SpacialSupport ) ;  
        
	void  _getUpdate( double * object, double * last_object, unsigned char * SpacialSupport ) ;
	void  _getUpdate( float  * object, float  * last_object, unsigned char * SpacialSupport ) ;     
} ;


#endif    /*   #include "deconvolver.h"   */
