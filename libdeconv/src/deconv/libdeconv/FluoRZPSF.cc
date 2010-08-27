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
 * Filename:  FluoRZPSF.cc
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "CSlice.h"
#include "FluoRZPSF.h" 


FluoRZPSF::FluoRZPSF( double NA, double WL, double RI, int DimR, int Sections, double RadialCalibration, double SectioningConstant )
{
	init( NA, WL, RI, DimR, Sections, RadialCalibration, SectioningConstant ) ;
}



void FluoRZPSF::init( double NA, double WL, double RI, int DimR, int Sections, double RadialCalibration, double SectioningConstant )
{	
	_init( NA, WL, RI ) ;
	
	if( DimR >= Samples_LowerLimit ) 
	{
		_DimR = DimR ;
	}
	else	throw SamplesError( DimR ) ;
	
	if( Sections >= Samples_LowerLimit ) 
	{
		_Sections = Sections ;
	}
	else	throw SamplesError( Sections ) ;
	
	if( SectioningConstant >= Sectioning_LowerLimit && SectioningConstant <= Sectioning_UpperLimit ) 
	{
		_DZ = SectioningConstant ;
	}
	else	throw SectioningError( SectioningConstant ) ;
	   
	if( RadialCalibration >= Calibration_LowerLimit && RadialCalibration <= Calibration_UpperLimit ) 
	{
		_DR = RadialCalibration ;
	}                        
	else	throw CalibrationError( RadialCalibration ) ;
}



int FluoRZPSF::maxDimensionX( double CalibrationX ) 
{ 
	return ( 2 * ( (int) floor( (double)_DimR * _DR / CalibrationX / 1.5 ) - 1 ) ) ; 
}



int FluoRZPSF::maxDimensionY( double CalibrationY ) 
{ 
	return ( 2 * ( (int) floor( (double)_DimR * _DR / CalibrationY / 1.5 ) - 1 ) ) ; 
}



int FluoRZPSF::maxSections( double SectioningConstant )
{
	return ( 2 * ( (int) floor( (double)_Sections * _DZ / SectioningConstant ) - 1 ) ) ;
}



double FluoRZPSF::minCalibration( int Points2Sum )  
{ 
	return ( _DR * ((double) Points2Sum) ) ; 
}



void FluoRZPSF::read( const char * filehead )
{
	char   filename[256] ;
	char   line[1024] ;
	double params[15] ;
	FILE * fp ;
	
	sprintf( filename, "%s.rzh", filehead ) ;
	fp = fopen( filename, "rt" ) ;	
	if( fp )
	{
		int i = 0 ;
		while( !feof( fp )  )
		{
			fgets( line, 1024, fp ) ;
			std::string s = std::string( line ) ;
			if( s.size() > 10 && s[0] != '#' ) 
			{
				sscanf( line, "%lf", &params[i] ) ;
				i++ ;
			} 
		}
		fclose( fp ) ;
		if( i >= 9 )
		{
			init( params[6], params[7], params[8], (int) params[2], (int) params[3], params[0], params[1] ) ;
			if( fabs( _NyqDXY - params[4] ) > 0.0001 || fabs( _NyqDF - params[5] ) > 0.0001 ) throw RZHError2( filename ) ;
			if( i > 9  ) setCoverSlip_Mismatch( params[12], params[11], params[10], params[9] ) ;
			if( i > 13 ) setImmersionMedium_Mismatch( params[13], params[8], params[14] ) ;
		}
		else	throw RZHError1( filename ) ;			
	}
	else	throw ErrnoError( std::string(filename) ) ;
	
	sprintf( filename, "%s.rzd", filehead ) ;
	fp = fopen( filename, "rb" ) ;
	if( fp )
	{
		DoubleArray temp( new double[ _DimR * _Sections ] ) ;
		_RZpsf = temp ;
		
		if( (int) fread( _RZpsf.get(), sizeof(double), (_DimR * _Sections), fp ) != (_DimR * _Sections) )
		{
			throw ReadDataError( std::string(filename) ) ;	
		}
		fclose( fp ) ;
	}
	else	throw ErrnoError( std::string(filename) ) ;
}



void FluoRZPSF::save( const char * filehead )
{
	char   filename[256] ;
	FILE * fp ;
	
	if( _RZpsf )
	{
		sprintf( filename, "%s.rzh", filehead ) ;
		fp = fopen( filename, "wt" ) ;
		if( fp )
		{			
			fprintf( fp, "%10.4f -> Radial Calibration Constant (um) in the object space.\n", _DR ) ;
			fprintf( fp, "%10.4f -> Optial Sectioning  Constant (um) along the optical axis.\n", _DZ ) ;		
			fprintf( fp, "%10d -> Radial   Samples.\n", _DimR ) ;
			fprintf( fp, "%10d -> Optical Sections.\n", _Sections  ) ;
			fprintf( fp, "\n" ) ;
			
			_exportCommon( fp ) ;
			fclose( fp ) ;
		}		
		else	throw ErrnoError( std::string(filename) ) ;
	
		sprintf( filename, "%s.rzd", filehead ) ;
		fp = fopen( filename, "wb" ) ;
		if( fp )
		{
			if( (int) fwrite( _RZpsf.get(), sizeof(double), (_DimR * _Sections), fp ) != (_DimR * _Sections) )
			{
				throw WriteDataError( std::string(filename) ) ;	
			}
			fclose( fp ) ;
		}
		else	throw ErrnoError( std::string(filename) ) ;
	}
	else	throw PSFError( 3 ) ;
}



void FluoRZPSF::create( bool Check )
{
	double defocus, r ;
	time_t t0, t1 ;
	
	if( _NA > 0.0 )
	{
		double NNA = _NA * _NA ;
		
		DoubleArray temp( new double[ _DimR * _Sections ] ) ;
		_RZpsf = temp ;
		
		gsl_integration_workspace * ws = gsl_integration_workspace_alloc( IntegrationLimit ) ;
	
		if( Check ) std::cout << " FluoRZPSF::create starts generating a RZ PSF : " << _DimR << "x" << _Sections << " ... \n" ;
	
		for( int j = 0 ; j < _Sections ; j++ )
		{
			if( Check ) time( &t0 ) ;
			defocus = ((double) j) * _DZ ;
			for( int i = 0 ; i < _DimR ; i++ )
			{
				r = ((double) i) * _DR ;
				_RZpsf[ i + j * _DimR ] = _IntegralPSF( ws, 0.0, NNA, r, defocus ) ;
			}
			if( Check ) 
			{
				time( &t1 ) ;
				std::cout << " Plane " << j << " processed and took " << difftime(t1, t0) << " seconds.\n" ;
			}
		}
        
		if( Check ) std::cout << " FluoRZPSF::create() completes generating the RZ_PSF.\n" ;
	}
	else	throw PSFError( 2 ) ;      
}



int FluoRZPSF::get3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, float * psf, const char * file, int Points2Sum )
{
	int    half_nx, half_ny, half_nz, ndz, sum_nx, sum_ny, half_sum_nx, half_sum_ny, half_Points2Sum, iii, jjj ;
	double dxx, dyy, dzz ; 
	double ddxx, ddyy, r, weight ;
	
	iii = _check3Dpsf( nx, ny, nz, dx, dy, dz, Points2Sum ) ;
	
	if( iii > 0 ) return iii ;
	else
	{	
		dxx = dx / ((double) Points2Sum ) ;
		dyy = dy / ((double) Points2Sum ) ;
		dzz = dz / _DZ ;
	
		DoubleArray x( new double[ _DimR ] ) ;
		DoubleArray y( new double[ _DimR ] ) ;
		
		for( int i = 0 ; i < _DimR ; i++ ) x[i] = ((double) i) * _DR ;
		
		ddxx = dxx * dxx ;
		ddyy = dyy * dyy ;
		half_nx = nx / 2 ;
		half_ny = ny / 2 ;
		half_nz = nz / 2 ;
		sum_nx  = nx * Points2Sum ;
		sum_ny  = ny * Points2Sum ;
		half_sum_nx = sum_nx / 2 ;
		half_sum_ny = sum_ny / 2 ;
		half_Points2Sum = ( Points2Sum - 1 ) / 2 ;
		ndz  = (int) (floor( dzz + 0.5 )) ;
	
		SingleSlice sum_slice( sum_nx, sum_ny ) ;
		SingleSlice slice( nx, ny ) ;
	
		gsl_interp_accel * acc = gsl_interp_accel_alloc() ;
		gsl_spline * spline = gsl_spline_alloc( gsl_interp_cspline, _DimR ) ;
	 
		for( int k = 0 ; k <= half_nz ; k++ )
		{
			slice.fillslice( 0.0 ) ;
			for( int i = 0 ; i < _DimR ; i++ ) y[i] = _RZpsf[ i + k * ndz * _DimR ] ;
			gsl_spline_init( spline, x.get(), y.get(), _DimR ) ;
			for( int j = 0 ; j <= half_sum_ny ; j++ )
				for( int i = 0 ; i <= half_sum_nx ; i++ )
				{
					r = sqrt( ((double)(i*i)) * ddxx + ((double)(j*j)) * ddyy ) ;
					weight = gsl_spline_eval( spline, r, acc ) ;
					if( j < half_sum_ny && i < half_sum_nx ) sum_slice( i, j )               = weight ;
					if( i > 0 )                 		 sum_slice( sum_nx-i, j )        = weight ;
					if( j > 0 )                  		 sum_slice( i, sum_ny-j )        = weight ;
					if( i > 0 && j > 0 )         		 sum_slice( sum_nx-i, sum_ny-j ) = weight ;
				}
			for( int j = 0 ; j < ny ; j++ )
				for( int i = 0 ; i < nx ; i++ )
				{
					for( int jj = Points2Sum*j-half_Points2Sum ; jj <= Points2Sum*j+half_Points2Sum ; jj++ )
					{
						if( jj < 0 )          	jjj = -jj ;
						else if( jj >= sum_ny ) jjj = 2 * sum_ny - jj - 2 ;
						else                    jjj = jj ;
						for( int ii = Points2Sum*i-half_Points2Sum ; ii <= Points2Sum*i+half_Points2Sum ; ii++ )
						{
							if( ii < 0 )          	iii = -ii ;
							else if( ii >= sum_nx ) iii = 2 * sum_ny - ii - 2 ;
							else            	iii = ii ;
							slice( i, j ) += sum_slice( iii, jjj ) ;
						}
					}
				}
			if( k == 0 || k == half_nz ) 
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (nz-k) * slice.size() ] = (slice.data())[i] ;
				}
			}
		}
        
		weight = 0.0 ;
		for( int i = 0 ; i < nx * ny * nz ; i++ ) weight += psf[i] ;
		for( int i = 0 ; i < nx * ny * nz ; i++ ) psf[i] /= weight ;
		
		if( file )
		{
			FILE* fp = fopen( file, "wt" ) ;
			if( fp )
			{
				_exportCommon( fp ) ;	
						
				fprintf( fp, "%10.4f -> Radial Calibration Constant (um) in the object space.\n", _DR ) ;
				fprintf( fp, "%10.4f -> Optial Sectioning  Constant (um) along the optical axis.\n", _DZ ) ;		
				fprintf( fp, "%10d -> Radial   Samples.\n", _DimR ) ;
				fprintf( fp, "%10d -> Optical Sections.\n", _Sections  ) ;
				fprintf( fp, "\n" ) ;
				
				fprintf( fp, "%10.4f -> Calibration Constant (um) along the X-axis in the object space.\n", dx ) ;
				fprintf( fp, "%10.4f -> Calibration Constant (um) along the Y-axis in the object space.\n", dy ) ;
				fprintf( fp, "%10.4f -> Sectioning  Constant (um) along the optical axis.\n", dz ) ;	
				fprintf( fp, "\n" ) ;
		
				fprintf( fp, "%10d -> Dimension of the PSF along the X-axis.\n", nx  ) ;
				fprintf( fp, "%10d -> Dimension of the PSF along the Y-axis.\n", ny  ) ;
				fprintf( fp, "%10d -> Sections  of the PSF along the optical axis.\n", nz  ) ;
				fprintf( fp, "\n" ) ;
				fclose( fp ) ;				
			}
		}

		gsl_spline_free( spline ) ;
		gsl_interp_accel_free( acc ) ;
		
		return 0 ;
	}
}



int FluoRZPSF::get3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, double * psf, const char * file, int Points2Sum )
{
	int    half_nx, half_ny, half_nz, ndz, sum_nx, sum_ny, half_sum_nx, half_sum_ny, half_Points2Sum, iii, jjj ;
	double dxx, dyy, dzz ;
	double ddxx, ddyy, r, weight ;
	
	iii = _check3Dpsf( nx, ny, nz, dx, dy, dz, Points2Sum ) ;
		
	if( iii > 0 ) return iii ;
	else
	{	
		dxx = dx / ((double) Points2Sum ) ;
		dyy = dy / ((double) Points2Sum ) ;
		dzz = dz / _DZ ;
	
		DoubleArray x( new double[ _DimR ] ) ;
		DoubleArray y( new double[ _DimR ] ) ;

		for( int i = 0 ; i < _DimR ; i++ ) x[i] = ((double) i) * _DR ;
		
		ddxx = dxx * dxx ;
		ddyy = dyy * dyy ;
		half_nx = nx / 2 ;
		half_ny = ny / 2 ;
		half_nz = nz / 2 ;
		sum_nx  = nx * Points2Sum ;
		sum_ny  = ny * Points2Sum ;
		half_sum_nx = sum_nx / 2 ;
		half_sum_ny = sum_ny / 2 ;
		half_Points2Sum = ( Points2Sum - 1 ) / 2 ;
		ndz  = (int) (floor( dzz + 0.5 )) ;	
	
		DoubleSlice sum_slice( sum_nx, sum_ny ) ;
		DoubleSlice slice( nx, ny ) ;
	
		gsl_interp_accel * acc = gsl_interp_accel_alloc() ;
		gsl_spline * spline = gsl_spline_alloc( gsl_interp_cspline, _DimR ) ;
	 
		for( int k = 0 ; k <= half_nz ; k++ )
		{
			slice.fillslice( 0.0 ) ;
			for( int i = 0 ; i < _DimR ; i++ ) y[i] = _RZpsf[ i + k * ndz * _DimR ] ;
			gsl_spline_init( spline, x.get(), y.get(), _DimR ) ;
			for( int j = 0 ; j <= half_sum_ny ; j++ )
				for( int i = 0 ; i <= half_sum_nx ; i++ )
				{
					r = sqrt( ((double)(i*i)) * ddxx + ((double)(j*j)) * ddyy ) ;
					weight = gsl_spline_eval( spline, r, acc ) ;
					if( j < half_sum_ny && i < half_sum_nx ) sum_slice( i, j )               = weight ;
					if( i > 0 )                 		 sum_slice( sum_nx-i, j )        = weight ;
					if( j > 0 )                  		 sum_slice( i, sum_ny-j )        = weight ;
					if( i > 0 && j > 0 )         		 sum_slice( sum_nx-i, sum_ny-j ) = weight ;
				}
			for( int j = 0 ; j < ny ; j++ )
				for( int i = 0 ; i < nx ; i++ )
				{
					for( int jj = Points2Sum*j-half_Points2Sum ; jj <= Points2Sum*j+half_Points2Sum ; jj++ )
					{
						if( jj < 0 )          	jjj = -jj ;
						else if( jj >= sum_ny ) jjj = 2 * sum_ny - jj - 2 ;
						else                    jjj = jj ;
						for( int ii = Points2Sum*i-half_Points2Sum ; ii <= Points2Sum*i+half_Points2Sum ; ii++ )
						{
							if( ii < 0 )          	iii = -ii ;
							else if( ii >= sum_nx ) iii = 2 * sum_ny - ii - 2 ;
							else            	iii = ii ;
							slice( i, j ) += sum_slice( iii, jjj ) ;
						}
					}
				}
			if( k == 0 || k == half_nz ) 
			{
				for( int i = 0 ; i < slice.size() ; i++ ) psf[ i + k * slice.size() ] = (slice.data())[i] ;
			}
			else
			{
				for( int i = 0 ; i < slice.size() ; i++ ) 
				{
					psf[ i + k * slice.size() ] = (slice.data())[i] ; 
					psf[ i + (nz-k) * slice.size() ] = (slice.data())[i] ;
				}
			}
		}
        
		weight = 0.0 ;
		for( int i = 0 ; i < nx * ny * nz ; i++ ) weight += psf[i] ;
		for( int i = 0 ; i < nx * ny * nz ; i++ ) psf[i] /= weight ;
		
		if( file )
		{
			FILE* fp = fopen( file, "wt" ) ;
			if( fp )
			{
				_exportCommon( fp ) ;	
						
				fprintf( fp, "%10.4f -> Radial Calibration Constant (um) in the object space.\n", _DR ) ;
				fprintf( fp, "%10.4f -> Optial Sectioning  Constant (um) along the optical axis.\n", _DZ ) ;		
				fprintf( fp, "%10d -> Radial   Samples.\n", _DimR ) ;
				fprintf( fp, "%10d -> Optical Sections.\n", _Sections  ) ;
				fprintf( fp, "\n" ) ;
				
				fprintf( fp, "%10.4f -> Calibration Constant (um) along the X-axis in the object space.\n", dx ) ;
				fprintf( fp, "%10.4f -> Calibration Constant (um) along the Y-axis in the object space.\n", dy ) ;
				fprintf( fp, "%10.4f -> Sectioning  Constant (um) along the optical axis.\n", dz ) ;	
				fprintf( fp, "\n" ) ;
		
				fprintf( fp, "%10d -> Dimension of the PSF along the X-axis.\n", nx ) ;
				fprintf( fp, "%10d -> Dimension of the PSF along the Y-axis.\n", ny  ) ;
				fprintf( fp, "%10d -> Sections  of the PSF along the optical axis.\n", nz  ) ;
				fprintf( fp, "\n" ) ;
				fclose( fp ) ;				
			}
		}
        
		gsl_spline_free( spline ) ;
		gsl_interp_accel_free( acc ) ;

		return 0 ;
	}
}



/* private functions */

int FluoRZPSF::_check3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, int Points2Sum )
{
	if( _RZpsf )
	{
		if( nx%2 != 0 || nx < Samples_LowerLimit || nx > maxDimensionX( dx ) )
		{
			std::cout << " Failed to get 3Dpsf for its X-dimension must be even and lager than "
			          << Samples_LowerLimit << " and smaller than " << maxDimensionX( dx ) << ".\n" ;
			return 1 ;
		}
	
		if( ny%2 != 0 || ny < Samples_LowerLimit || ny > maxDimensionY( dy ) )
		{
			std::cout << " Failed to get 3Dpsf for its Y-dimension must be even and lager than "
			          << Samples_LowerLimit << " and smaller than " << maxDimensionX( dy ) << ".\n" ;
			return 2 ;		
		}
	    
		if( nz%2 != 0 || nz < Samples_LowerLimit || nz > maxSections( dz ) )
		{
			std::cout << " Failed to get 3Dpsf for its Z-dimension must be even and lager than "
			          << Samples_LowerLimit << " and smaller than " << maxSections( dz ) << ".\n" ;
			return 3 ;
		}	
	
		if( Points2Sum%2 == 0 || Points2Sum < 1 ) 
		{
			std::cout << " Failed to get 3Dpsf for the compensation parameter <Points2Sum> must be an odd integer.\n" ;
			return 4 ;
		}
	
		if( dx < minCalibration( Points2Sum ) || dy < minCalibration( Points2Sum ) ) 
		{
			std::cout << " Failed to get 3Dpsf for its calibrations are so small that RZpsf can not generate it.\n" 
			          << " Try to lower the compensation parameter <Points2Sum>.\n" ;
			return 5 ;
		}
	
		if( ( dz/_DZ - floor(dz/_DZ + 0.5) ) > DiffEpsilon || dz < _DZ )
		{
			std::cout << " Failed to get 3Dpsf for its sectioning constant must be N*" 
			          << _DZ << " where N is an positive integer.\n" ;
			return 6 ;
		}
	}
	else	throw PSFError( 3 ) ;
	
	return 0 ;
}
