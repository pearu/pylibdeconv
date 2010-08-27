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
 * Filename:  LWCGdeconvolver.h
 */
 

#ifndef LWCGDECONVOLVER_H
#define LWCGDECONVOLVER_H


#include "deconvolver.h"
#include "FFTW3fft.h"


/*
 *	==============================================================================
 *	LWCGdeconvolver is a base class for the LWdeconvoler and CGdeconvolver classes 
 *	since a same pre-conditioning technique is required for both of them.
 *	==============================================================================
 *
 *	------------------------------------------------------------------------------------------
 *	Pre-Conditioning: <_ConditioningIteration>, <_ConditioningTolerance>, <_ConditioningValue>
 *	-------------------------------------------------------------------------------------------
 *		
 *	<_ConditioningIteration> : it is the number of iterations used in golden search for <_ConditioningValue>;
 *	                           the larger <_ConditioningIteration>, the finer <_ConditioningValue> obtained,
 *	                           but it will take more time; its default value is 5. 
 *
 *	<_ConditioningTolerance> : it is the tolerance used in golden search for <_ConditioningValue>; 
 *	                           the smaller <_ConditioningTolerance>, the finer <_ConditioningValue> obtained,
 *	                           but it will take more time; its default value is 0.1. 
 *
 *	<_ConditioningValue>     : it is the conditioning value used for the pre-conditioning; its default value is 1.0e-8.
 *		
 *	If <_ConditioingIteration> is not set to be zero by the user, 
 *	<_ConditioningValue> used for the pre-conditioning will be calculated in {LW/CG}deconvolver::run(). 
 *	
 *	If <_ConditioingIteration> is set to be zero by the user, 
 *	<_ConditioningValue> will not be calculated be calculated in {LW/CG}deconvolver::run(), 
 *	instead {LW/CG}deconvolver::run() will use <_ConditioningValue> set by user to do pre-conditioning; 
 *	Make sure that you can set <_ConditioningValue> correctly.
 */


#define ConditioningIterationLimit 50
#define ConditioningValueLimit     1.0E-10
#define ConditioningToleranceLimit 1.0E-6


class ConditioningIterationError : public Error
{
        public:
        ConditioningIterationError( unsigned char iter )
        {
                _error << " Conditioning_Iteration Setup Error ( it must be no larger than " 
                       << ConditioningIterationLimit << " ) :\n"
                       << " Conditioning_Iteration was set -> " << iter << "\n" ;
        }
} ;


class ConditioningValueError : public Error
{
        public:
        ConditioningValueError( double cv )
        {
                _error << " Conditioning_Value Setup Error ( it must be larger than " 
                       << ConditioningValueLimit << " ) :\n"
                       << " Conditioning_Value was set -> " << cv << "\n" ;
        }
} ;


class ConditioningToleranceError : public Error
{
        public:
        ConditioningToleranceError( double ct )
        {
                _error << " Conditioning_Tolerance Setup Error ( it must be larger than " 
                       << ConditioningToleranceLimit << " ) :\n"
                       << " Conditioning_Tolerance was set -> " << ct << "\n" ;
        }
} ;

               
class LWCGdeconvolver : public deconvolver
{
	public:	
	virtual         ~LWCGdeconvolver() {}
	
	
	/*
	 *	Get protected members
	 *	ConditioningIteration() returns <_ConditioningIteration> described above.
	 *	ConditioningTolerance() returns <_ConditioningTolerance> described above.
	 *	ConditioningValue()     returns <_ConditioningValue>     described above.
	 */
	unsigned int  ConditioningIteration()  { return _ConditioningIteration ; }
	double        ConditioningTolerance()  { return _ConditioningTolerance ; }
	double        ConditioningValue()      { return _ConditioningValue ;     }
	
	
	/*
	 *	Set <_ConditioningIteration> described above
	 *	Input:
	 *		iter, it is the number of iterations and its default value is 5.
	 *	Throw:
	 *		throw an error if the input is out of the pre-defiend range.
	 */
	void    setConditioningIteration( unsigned int iter = 5 ) ;
	
	
	/*
	 *	Set <_ConditioningTolerance> described above
	 *	Input:
	 *		ct, it is the tolerance and its default value is 0.1.
	 *	Throw:
	 *		throw an error if the input is out of the pre-defined range.
	 */
	void    setConditioningTolerance( double ct = 0.1 ) ;
	
	
	/*
	 *	Set <_ConditioningValue> described above
	 *	Input:
	 *		cv, it is the conditioning value and its default value is 1.0e-8.
	 *	Throw:
	 *		throw an error if the input is out of the pre-defined range.
	 */
	void    setConditioningValue( double cv = 1.0e-8 ) ;
        

	protected:        
	unsigned int    _ConditioningIteration ;
	double          _ConditioningValue ;
	double          _ConditioningTolerance ;	
	FFTW3_FFT::Ptr  _FFTplanf ;
	FFTW3_FFT::Ptr  _FFTplanb ;
	
	void   	_printConditioning( double cv, double likelihood ) ;
	
	void    _exportLWCG( FILE * fp ) ;
	
	double  _getLikelihood( int size, double * image_re, double * image_im,  double * psf_re,
	                        double * psf_im,   double * object_re, double * object_im ) ;
                               
	double  _getLikelihood( int size, float  * image_re, float  * image_im,  float  * psf_re,
	                        float  * psf_im,   float  * object_re, float  * object_im ) ;
                                                  
	void    _update( int size, double cv, double * object_re, double * object_im, 
	                 double * image_re, double * image_im, double * psf_re, double * psf_im, double * otf ) ;
                                   
	void    _update( int size, double cv, float * object_re, float * object_im, 
	                 float * image_re, float * image_im, float * psf_re, float * psf_im, float * otf ) ;
			
	double  _getSumLikelihood( int size, double cv, double * object0, double * object, 
	                           double * object_re, double * object_im, double * image_re, double * image_im, 
	                           double * psf_re, double * psf_im, double * otf, unsigned char * SpacialSupport ) ;

	double  _getSumLikelihood( int size, double cv, float * object0, float * object, 
	                           float * object_re, float * object_im, float * image_re, float * image_im, 
	                           float * psf_re, float * psf_im, float * otf, unsigned char * SpacialSupport ) ;

	double  _runConditioning( int size, double * object0, double * object, 
	                          double * object_re, double * object_im, double * image_re, double * image_im, 
	                          double * psf_re, double * psf_im, double * otf, unsigned char * SpacialSupport ) ;
				 
	double  _runConditioning( int size, float * object0, float * object, 
	                          float * object_re, float * object_im, float * image_re, float * image_im, 
	                          float * psf_re, float * psf_im, float * otf, unsigned char * SpacialSupport ) ;
} ;


#endif   /*   #include "LWCGdeconvolver.h"   */
