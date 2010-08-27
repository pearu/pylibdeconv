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
 * Filename:  Fluo3DPSF.cc
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CSlice.h"
#include "Fluo3DPSF.h" 



Fluo3DPSF::Fluo3DPSF( double NA, double WL, double RI, double CalibrationX, double CalibrationY, double SectioningConstant )
{
	init( NA, WL, RI, CalibrationX, CalibrationY, SectioningConstant ) ;
}



void Fluo3DPSF::init( double NA, double WL, double RI, double CalibrationX, double CalibrationY, double SectioningConstant )
{
	_init( NA, WL, RI ) ;
	
	if( CalibrationX >= Calibration_LowerLimit && CalibrationX <= Calibration_UpperLimit ) 
	{
		_DX = CalibrationX ;     
	}                   
	else	throw CalibrationError( CalibrationX ) ;
		
	if( CalibrationY >= Calibration_LowerLimit && CalibrationY <= Calibration_UpperLimit ) 
	{
		_DY = CalibrationY ;
	}                        
	else	throw CalibrationError( CalibrationY ) ;
		
	if( SectioningConstant >= Sectioning_LowerLimit && SectioningConstant <= Sectioning_UpperLimit ) 
	{
		_DZ = SectioningConstant ;
	}
	else	throw SectioningError( SectioningConstant ) ;

	_DimX = 0 ;
	_DimY = 0 ;
	_DimZ = 0 ; 
}



void Fluo3DPSF::exportProfile( const char * filename )
{
	if( _NA > 0.0 )
	{
		FILE * fp = fopen( filename, "wt" ) ;
		if( fp )
		{	
			_exportCommon( fp ) ;
			
			fprintf( fp, "%10.4f -> Calibration Constant (um) along the X-axis in the object space.\n", _DX ) ;
			fprintf( fp, "%10.4f -> Calibration Constant (um) along the Y-axis in the object space.\n", _DY ) ;
			fprintf( fp, "%10.4f -> Sectioning  Constant (um) along the optical axis.\n", _DZ ) ;	
			fprintf( fp, "\n" ) ;
		
			fprintf( fp, "%10d -> Dimension of the PSF along the X-axis.\n", _DimX  ) ;
			fprintf( fp, "%10d -> Dimension of the PSF along the Y-axis.\n", _DimY  ) ;
			fprintf( fp, "%10d -> Sections  of the PSF along the optical axis.\n", _DimZ  ) ;
			fprintf( fp, "\n" ) ;
			
			fclose( fp );
		}
		else	throw ErrnoError( std::string(filename) ) ;
	}
	else	throw PSFError( 1 ) ;
}



void Fluo3DPSF::create( int DimX, int DimY, int DimZ, float * psf, bool Check )
{
	time_t t0, t1 ;
	double defocus, r, weight ;	
	double NNA = _NA * _NA ;
	double DDX = _DX * _DX ;
	double DDY = _DY * _DY ;
	int  HalfX = DimX / 2 ;
	int  HalfY = DimY / 2 ;
	int  HalfZ = DimZ / 2 ;

	_setDimensions( DimX, DimY, DimZ ) ;

	SingleSlice slice( _DimX, _DimY ) ;
	gsl_integration_workspace * ws = gsl_integration_workspace_alloc( IntegrationLimit ) ;
	
	if( Check ) std::cout << " Fluo3DPSF::create starts generating a 3D PSF : "
	                      << _DimX << "x" << _DimY << "x" << _DimZ << " ... \n" ;
	
	if( fabs( DDX - DDY ) < DiffEpsilon )
	{
		DoubleSlice val( HalfX+1, HalfY+1 ) ;
		for( int k = 0 ; k <= HalfZ ; k++ )
		{
			if( Check ) time( &t0 ) ;
			val.fillslice( -1.0 ) ;	
			defocus = ((double) k) * _DZ ;
			for( int j = 0 ; j <= HalfY ; j++ )
				for( int i = 0 ; i <= HalfX ; i++ )
				{
					if( j <= HalfX && i <= HalfY )
					{
						if( val( j, i ) < 0.0 )
						{
							r = sqrt( ((double)(i*i+j*j)) * DDX ) ;
							val( i, j ) = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
						}
						else	val( i, j ) = val( j, i ) ;
					}
					else
					{
						r = sqrt( ((double)(i*i+j*j)) * DDX ) ;
						val( i, j ) = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
					} 
					if( j < HalfY && i < HalfX ) slice( i, j )             = val( i, j ) ;
					if( i > 0 )                  slice( _DimX-i, j )       = val( i, j ) ;
					if( j > 0 )                  slice( i, _DimY-j )       = val( i, j ) ;
					if( i > 0 && j > 0 )         slice( _DimX-i, _DimY-j ) = val( i, j ) ;
				}
			if( k == 0 || k == HalfZ )
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Plane " << k << " processed and took " 
					          << difftime(t1, t0) << " seconds.\n" ;
				}
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (_DimZ-k) * slice.size() ] = (slice.data())[i] ;
				}
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Planes " << k << " and " << _DimZ-k << " processed and took " 
						  << difftime(t1, t0) << " seconds.\n" ;
				} 			
			}	
		}
	}
	else
	{
		for( int k = 0 ; k <= HalfZ ; k++ )
		{
			if( Check ) time( &t0 ) ;
			defocus = ((double) k) * _DZ ;
			for( int j = 0 ; j <= HalfY ; j++ )
				for( int i = 0 ; i <= HalfX ; i++ )
				{
					r = sqrt( ((double)(i*i)) * DDX + ((double)(j*j)) * DDY ) ;
					weight = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
					if( j < HalfY && i < HalfX ) slice( i, j )             = weight ;
					if( i > 0 )                  slice( _DimX-i, j )       = weight ;
					if( j > 0 )                  slice( i, _DimY-j )       = weight ;
					if( i > 0 && j > 0 )         slice( _DimX-i, _DimY-j ) = weight ;
				}
			if( k == 0 || k == HalfZ )
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Plane " << k << " processed and took " 
					          << difftime(t1, t0) << " seconds.\n" ;
		 		}
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (_DimZ-k) * slice.size() ] = (slice.data())[i] ;
				}
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Planes " << k << " and " << _DimZ-k << " processed and took " 
					          << difftime(t1, t0) << " seconds.\n" ;
				} 			
			}
		}
	}
	
	weight = 0.0 ;
	for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ ) weight += psf[i] ;
	for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ ) psf[i] /= weight ;
	
	if( Check ) std::cout << " Fluo3DPSF::create() completes generating the 3D PSF.\n" ;
}



void Fluo3DPSF::create( int DimX, int DimY, int DimZ, double * psf, bool Check )
{
	time_t t0, t1 ;	
	double defocus, r, weight ;	
	double NNA = _NA * _NA ;
	double DDX = _DX * _DX ;
	double DDY = _DY * _DY ;
	int  HalfX = DimX / 2 ;
	int  HalfY = DimY / 2 ;
	int  HalfZ = DimZ / 2 ;
	
	_setDimensions( DimX, DimY, DimZ ) ;
	
	DoubleSlice slice( _DimX, _DimY ) ;
	gsl_integration_workspace * ws = gsl_integration_workspace_alloc( IntegrationLimit ) ;
	
	if( Check ) std::cout << " Fluo3DPSF::create starts generating a 3D PSF : "
	                      << _DimX << "x" << _DimY << "x" << _DimZ << " ... \n" ;
	
	if( fabs( DDX - DDY ) < DiffEpsilon )
	{
		DoubleSlice val( HalfX+1, HalfY+1 ) ;
		for( int k = 0 ; k <= HalfZ ; k++ )
		{
			if( Check ) time( &t0 ) ;
			val.fillslice( -1.0 ) ;	
			defocus = ((double) k) * _DZ ;
			for( int j = 0 ; j <= HalfY ; j++ )
				for( int i = 0 ; i <= HalfX ; i++ )
				{
					if( j <= HalfX && i <= HalfY )
					{
						if( val( j, i ) < 0.0 )
						{
							r = sqrt( ((double)(i*i+j*j)) * DDX ) ;
							val( i, j ) = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
						}
						else	val( i, j ) = val( j, i ) ;
					}
					else
					{
						r = sqrt( ((double)(i*i+j*j)) * DDX ) ;
						val( i, j ) = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
					} 
					if( j < HalfY && i < HalfX ) slice( i, j )             = val( i, j ) ;
					if( i > 0 )                  slice( _DimX-i, j )       = val( i, j ) ;
					if( j > 0 )                  slice( i, _DimY-j )       = val( i, j ) ;
					if( i > 0 && j > 0 )         slice( _DimX-i, _DimY-j ) = val( i, j ) ;
				}
			if( k == 0 || k == HalfZ )
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Plane " << k << " processed and took " 
						  << difftime(t1, t0) << " seconds.\n" ;
				}
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (_DimZ-k) * slice.size() ] = (slice.data())[i] ;
				}
				if( Check ) 
				{
					time( &t1 ) ;
					std::cout << " Planes " << k << " and " << _DimZ-k << " processed and took " 
					          << difftime(t1, t0) << " seconds.\n" ;			
				}
			}	
		}
	}
	else
	{
		for( int k = 0 ; k <= HalfZ ; k++ )
		{
			if( Check ) time( &t0 ) ;
			defocus = ((double) k) * _DZ ;
			for( int j = 0 ; j <= HalfY ; j++ )
				for( int i = 0 ; i <= HalfX ; i++ )
				{
					r = sqrt( ((double)(i*i)) * DDX + ((double)(j*j)) * DDY ) ;
					weight = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
					if( j < HalfY && i < HalfX ) slice( i, j )             = weight ;
					if( i > 0 )                  slice( _DimX-i, j )       = weight ;
					if( j > 0 )                  slice( i, _DimY-j )       = weight ;
					if( i > 0 && j > 0 )         slice( _DimX-i, _DimY-j ) = weight ;
				}
			if( k == 0 || k == HalfZ )
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
				if( Check )
				{
					time( &t1 ) ;
					std::cout << " Plane " << k << " processed and took " 
						  << difftime(t1, t0) << " seconds.\n" ;	
				}
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (_DimZ-k) * slice.size() ] = (slice.data())[i] ;
				}
				if( Check )
				{
					time( &t1 ) ;
					std::cout << " Planes " << k << " and " << _DimZ-k << " processed and took " 
					          << difftime(t1, t0) << " seconds.\n" ;
				}	 			
			}
		}
	}
	
	weight = 0.0 ;
	for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ ) weight += psf[i] ;
	for( int i = 0 ; i < _DimX * _DimY * _DimZ ; i++ ) psf[i] /= weight ;
	
	if( Check ) std::cout << " Fluo3DPSF::create() completes generating the 3D PSF.\n" ;
}



/* private function */

void Fluo3DPSF::_setDimensions( int DimX, int DimY, int DimZ )
{
	if( _NA > 0.0 )
	{
		if( (DimX%2) == 0 && DimX >= Samples_LowerLimit )
		{
			_DimX = DimX ;
		}
		else	throw SamplesError( DimX ) ;
	
		if( (DimY%2) == 0 && DimY >= Samples_LowerLimit ) 
		{
			_DimY = DimY ;
		}
		else	throw SamplesError( DimY ) ;
	
		if( (DimZ%2) == 0 && DimZ >= Samples_LowerLimit ) 
		{
			_DimZ = DimZ ;
		}
		else	throw SamplesError( DimZ ) ;
	} 
	else	throw PSFError( 1 ) ;
}
