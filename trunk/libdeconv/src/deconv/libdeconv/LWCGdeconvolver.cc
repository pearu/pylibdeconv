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
 * Filename:  LWCGdeconvolver.cc
 */
 

#include "math.h"
#include "LWCGdeconvolver.h"


/* public functions */ 

void LWCGdeconvolver::setConditioningIteration( unsigned int iter )
{ 
	if( iter <= ConditioningIterationLimit ) _ConditioningIteration = iter ; 
	else                                     throw ConditioningIterationError( iter ) ;
}


        
void LWCGdeconvolver::setConditioningTolerance( double ct )
{
	if( ct >= ConditioningToleranceLimit ) _ConditioningTolerance = ct ;
	else                                   throw ConditioningToleranceError( ct ) ;
}



void LWCGdeconvolver::setConditioningValue( double cv )
{
	if( cv >= ConditioningValueLimit ) _ConditioningValue = cv ;
	else                               throw ConditioningValueError( cv ) ;
}



/* protected functions */

void LWCGdeconvolver::_printConditioning( double cv, double likelihood ) 
{
	printf( " --> Conditioning value = %9.6f -> Likelihood = %12.6e\n", cv, likelihood ) ; 
}



void LWCGdeconvolver::_exportLWCG( FILE * fp )
{
	_exportCommon( fp ) ;
	fprintf( fp, "%d -> Conditioning_Iteration of each conditioning process before the deconvolution loop.\n", _ConditioningIteration ) ;
	fprintf( fp, "%e -> Conditioning_Tolerance in each conditioning process before the deconvolution loop.\n", _ConditioningTolerance ) ;
	fprintf( fp, "%e -> Conditioning_Value found after conditioning and used in the deconvolution loop\n", _ConditioningValue ) ;
	fprintf( fp, "\n" ) ;
}



double LWCGdeconvolver::_getLikelihood( int size, double * image_re, double * image_im, double * psf_re, 
                        double * psf_im, double * object_re, double * object_im )
{ 
	double likelihood = 0.0 ;
     	
	for( int i = 0 ; i < size ; i++ )
	{
		likelihood += ( ( image_re[i]-psf_re[i]*object_re[i]+psf_im[i]*object_im[i] ) *
		                ( image_re[i]-psf_re[i]*object_re[i]+psf_im[i]*object_im[i] ) +
		                ( image_im[i]-psf_re[i]*object_im[i]-psf_im[i]*object_re[i] ) *
		                ( image_im[i]-psf_re[i]*object_im[i]-psf_im[i]*object_re[i] ) ) ;
       		
	}
     	
	return likelihood ;
}



double LWCGdeconvolver::_getLikelihood( int size, float * image_re, float * image_im, float * psf_re, 
                        float * psf_im, float * object_re, float * object_im )
{ 
	double likelihood = 0.0 ;
	
	for( int i = 0 ; i < size ; i++ )
	{
		likelihood += ( ( image_re[i]-psf_re[i]*object_re[i]+psf_im[i]*object_im[i] ) *
		                ( image_re[i]-psf_re[i]*object_re[i]+psf_im[i]*object_im[i] ) +
		                ( image_im[i]-psf_re[i]*object_im[i]-psf_im[i]*object_re[i] ) *
		                ( image_im[i]-psf_re[i]*object_im[i]-psf_im[i]*object_re[i] ) ) ;
	}
     	
	return likelihood ;
}



void LWCGdeconvolver::_update( int size, double cv, double * object_re, double * object_im, 
                      double * image_re, double * image_im, double * psf_re, double * psf_im, double * otf )
{
	double temp1, temp2, temp3, temp4 ;
	
	for( int i = 0 ; i < size ; i++ )
	{
		       temp1 = otf[i] + cv ;
		       temp2 = ( image_re[i]*psf_re[i] + image_im[i]*psf_im[i] ) / temp1 ;
		       temp3 = ( image_im[i]*psf_re[i] - image_re[i]*psf_im[i] ) / temp1 ;
		       temp4 = otf[i] / temp1 ;
		object_re[i] = temp2 - temp4 * object_re[i] ;
		object_im[i] = temp3 - temp4 * object_im[i] ;
	}
	       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;
}



void LWCGdeconvolver::_update( int size, double cv, float * object_re, float * object_im, 
                      float * image_re, float * image_im, float * psf_re, float * psf_im, float * otf )
{
	float temp1, temp2, temp3, temp4 ;
	
	for( int i = 0 ; i < size ; i++ )
	{
		       temp1 = otf[i] + cv ;
		       temp2 = ( image_re[i]*psf_re[i] + image_im[i]*psf_im[i] ) / temp1 ;
		       temp3 = ( image_im[i]*psf_re[i] - image_re[i]*psf_im[i] ) / temp1 ;
		       temp4 = otf[i] / temp1 ;
		object_re[i] = temp2 - temp4 * object_re[i] ;
		object_im[i] = temp3 - temp4 * object_im[i] ;
	}
       	
	_FFTplanb->execute( object_re, object_im, object_re ) ;	
}



double LWCGdeconvolver::_getSumLikelihood( int size, double cv, double * object0, double * object, double * object_re, double * object_im,
                        double * image_re, double * image_im, double * psf_re, double * psf_im, double * otf, unsigned char * SpacialSupport )
{
	double sum_likelihood = 0.0 ;
	unsigned int iter = 0 ;

	for( int i = 0 ; i < _Space ; i++ ) object[i] = object0[i] ;
	_FFTplanf->execute( object, object_re, object_im ) ;

	while( iter < _ConditioningIteration )
	{
		_update( size, cv, object_re, object_im, image_re, image_im, psf_re, psf_im, otf ) ;

		for( int i = 0 ; i < _Space ; i++ )
		{
			object[i] += object_re[i] ;
			if( object[i] < 0.0 ) object[i] = 0.0 ;
		}
       		
		if( SpacialSupport != NULL )
		{
			for( int i = 0 ; i < _Space ; i++ ) object[i] *= ((double) SpacialSupport[i]) ;
		}
	
		_FFTplanf->execute( object, object_re, object_im ) ;

		sum_likelihood += _getLikelihood( size, image_re, image_im, psf_re, psf_im, object_re, object_im ) ;

		iter++ ;
	} 
    	   	
	return sum_likelihood ;
}



double LWCGdeconvolver::_getSumLikelihood( int size, double cv, float * object0, float * object, float * object_re, float * object_im,
                        float * image_re, float * image_im, float * psf_re, float * psf_im, float * otf, unsigned char * SpacialSupport )
{	
	unsigned int iter = 0 ;	
	double sum_likelihood = 0.0 ;

	for( int i = 0 ; i < _Space ; i++ ) object[i] = object0[i] ;
	_FFTplanf->execute( object, object_re, object_im ) ;

	while( iter < _ConditioningIteration )
	{
		_update( size, cv, object_re, object_im, image_re, image_im, psf_re, psf_im, otf ) ;
		
		for( int i = 0 ; i < _Space ; i++ )
		{
			object[i] += object_re[i] ;
			if( object[i] < 0.0 ) object[i] = 0.0 ;
		}
       		
		if( SpacialSupport != NULL )
		{
			for( int i = 0 ; i < _Space ; i++ ) object[i] *= ((double) SpacialSupport[i]) ;
		}

		_FFTplanf->execute( object, object_re, object_im ) ;

		sum_likelihood += _getLikelihood( size, image_re, image_im, psf_re, psf_im, object_re, object_im ) ;

		iter++ ;
	}
 	
	return sum_likelihood ;
}



double LWCGdeconvolver::_runConditioning( int size, double * object0, double * object, double * object_re, double * object_im, double * image_re,
                        double * image_im, double * psf_re, double * psf_im, double * otf, unsigned char * SpacialSupport )
{
	double R   = 0.61803399 ;
	double C   = 1.0 - R ;
	double likelihood = 1.0E+37 ;
	double last_likelihood = likelihood ;
	double x0, x1, x2, x3, f1, f2 ;

	x0 = 1.0 ;
	while( likelihood <= last_likelihood )
	{
		last_likelihood = likelihood ;
		x0 = x0 * 0.1 ;
		likelihood = _getSumLikelihood( size, x0, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x0, likelihood ) ;
	}        	       	       

	x1 = x0 * 10.0 ;
	x3 = x0 * 100.0 ;    	
	if( fabs(x3-x1) > fabs(x1-x0) )
	{
		f1 = last_likelihood ;
		x2 = x1 + C * ( x3 - x1 ) ;       		
		f2 = _getSumLikelihood( size, x2, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x2, f2 ) ;
	}
	else
	{
		x2 = x1 ;
		f2 = last_likelihood ;
		x1 = x2 - C * ( x2 - x0 ) ;
		f1 = _getSumLikelihood( size, x1, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x1, f1 ) ;
	}
       	
	while( fabs(x3-x0) > _ConditioningTolerance * ( fabs(x1) + fabs(x2) ) )
	{
		if ( f2 < f1 )
		{
			x0 = x1 ;
			x1 = x2 ;
			x2 = R * x1 + C * x3 ;
			f1 = f2 ;
			f2 = _getSumLikelihood( size, x2, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
			if( _CheckStatus ) _printConditioning( x2, f2 ) ;
		}
		else
		{ 
			x3 = x2 ; 
			x2 = x1 ;
			x1 = R * x2 + C * x0 ;
			f2 = f1 ;
			f1 = _getSumLikelihood( size, x1, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
			if( _CheckStatus ) _printConditioning( x1, f1 ) ;
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



double LWCGdeconvolver::_runConditioning( int size, float * object0, float * object, float * object_re, float * object_im, float * image_re, 
                        float * image_im, float * psf_re, float * psf_im, float * otf, unsigned char * SpacialSupport )
{
	double R   = 0.61803399 ;
	double C   = 1.0 - R ;
	double likelihood = 1.0E+37 ;
	double last_likelihood = likelihood ;
	double x0, x1, x2, x3, f1, f2 ;

	x0 = 1.0 ;
	while( likelihood <= last_likelihood )
	{
		last_likelihood = likelihood ;
		x0 = x0 * 0.1 ;
		likelihood = _getSumLikelihood( size, x0, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x0, likelihood ) ;
	}        	       	       

	x1 = x0 * 10.0 ;
	x3 = x0 * 100.0 ;    	
	if( fabs(x3-x1) > fabs(x1-x0) )
	{
		f1 = last_likelihood ;
		x2 = x1 + C * ( x3 - x1 ) ;       		
		f2 = _getSumLikelihood( size, x2, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x2, f2 ) ;
	}
	else
	{
		x2 = x1 ;
		f2 = last_likelihood ;
		x1 = x2 - C * ( x2 - x0 ) ;
		f1 = _getSumLikelihood( size, x1, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
		if( _CheckStatus ) _printConditioning( x1, f1 ) ;
	}
	   	
	while( fabs(x3-x0) > _ConditioningTolerance * ( fabs(x1) + fabs(x2) ) )
	{
		if( f2 < f1 )
		{
			x0 = x1 ;
			x1 = x2 ;
			x2 = R * x1 + C * x3 ;
			f1 = f2 ;
			f2 = _getSumLikelihood( size, x2, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
			if( _CheckStatus ) _printConditioning( x2, f2 ) ;
		}
		else
		{ 
			x3 = x2 ; 
			x2 = x1 ;
			x1 = R * x2 + C * x0 ;
			f2 = f1 ;
			f1 = _getSumLikelihood( size, x1, object0, object, object_re, object_im, image_re, image_im, psf_re, psf_im, otf, SpacialSupport ) ;
			if( _CheckStatus ) _printConditioning( x1, f1 ) ;
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

