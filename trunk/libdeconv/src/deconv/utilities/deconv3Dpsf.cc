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
 * Filename:  deconv3Dpsf.cc
 */
 

#include <stdlib.h>
#include "libdeconv/Fluo3DPSF.h"
#include "libdeconv/FluoRZPSF.h"
#include "libdeconv/CCube.h"


/*
 *	deconv3Dpsf is designed for creating a 3-D PSF using the Fluo3DPSF and CCube classes.
 *
 *	Eight input arguments are required by the program:
 *		<psf> <PSF.para> <DimensionX> <DimensionY> <DimensionZ> <CalibrationX> <CalibrationY> <SectioningZ>
 *
 *	The generated PSF will be saved using CCube::save():
 *		The header file : <psf>.hdr
 *		The data   file : <psf>.f32 if <DataType> is float. 
 *		                  <psf>.f64 if <DataType> is double.
 *
 *	The profile of the generated PSF will be saved in a text file : <psf>.log
 */


#define DataType           float                   
#define CheckPSFgeneration true


char usage[] = 
"$deconv3Dpsf <psf> <PSF.para> <DimensionX> <DimensionY> <DimensionZ> <CalibrationX> <CalibrationY> <SectioningZ>\n" ;


int main( int argc, char ** argv )
{
	Fluo3DPSF        psf ;
	CCube<DataType>  MyPsf ;
	char             filename[256] ;
	char             line[1024] ;
	double           Objective_Numerical_Aperture, Emission_Light_Wavelength, Act_Immersion_Medium_Refractive_Index ;
	double           Req_Immersion_Medium_Refractive_Index, Objective_Working_Distance ;
	double           Req_Cover_Refractive_Index, Act_Cover_Refractive_Index, Req_Cover_Thickness, Act_Cover_Thickness ;
	int              mode ;

	if( argc == 9 )
	{
		sprintf( filename, "%s.log", argv[1] ) ;
		
		FILE * fp = fopen( argv[2], "rt" ) ;
		if( fp ) 
		{
			fgets( line, 1024, fp ) ;
			sscanf( line, "%d", &mode ) ;
			if( mode >= 0 && mode <= 3 )
			{
				fgets( line, 1024, fp ) ;
				sscanf( line, "%lf", &Objective_Numerical_Aperture ) ;
				fgets( line, 1024, fp ) ;
				sscanf( line, "%lf", &Emission_Light_Wavelength ) ;
				fgets( line, 1024, fp ) ;
				sscanf( line, "%lf", &Act_Immersion_Medium_Refractive_Index ) ;	
				if( mode > 0 )
				{
					fgets( line, 1024, fp ) ;
					sscanf( line, "%lf", &Req_Immersion_Medium_Refractive_Index ) ;
					fgets( line, 1024, fp ) ;
					sscanf( line, "%lf", &Objective_Working_Distance ) ;
					if( mode > 1 )
					{
						fgets( line, 1024, fp ) ;
						sscanf( line, "%lf", &Act_Cover_Refractive_Index ) ;
						fgets( line, 1024, fp ) ;
						sscanf( line, "%lf", &Req_Cover_Refractive_Index ) ;
						fgets( line, 1024, fp ) ;
						sscanf( line, "%lf", &Act_Cover_Thickness ) ;
						fgets( line, 1024, fp ) ;
						sscanf( line, "%lf", &Req_Cover_Thickness ) ;
					}
				}								
			}
			else
			{
				std::cout << " program failed - check the PSF.para file.\n" ;
				return 1 ;				
			}
			fclose( fp ) ;
		}
		else
		{
			std::cout << " program failed to open the PSF.para file.\n" ;
			return 1 ;
		}
		
		int    dimx = atoi( argv[3] ) ;
		int    dimy = atoi( argv[4] ) ;
		int    dimz = atoi( argv[5] ) ;		
		double Calibration_Along_X_Dimension = atof( argv[6] ) ; 
		double Calibration_Along_Y_Dimension = atof( argv[7] ) ;
		double Optical_Sectioning_Constant   = atof( argv[8] ) ;
		
		MyPsf.init( dimx, dimy, dimz ) ;
		
		psf.init( Objective_Numerical_Aperture, Emission_Light_Wavelength, 
		          Act_Immersion_Medium_Refractive_Index, Calibration_Along_X_Dimension, 
		          Calibration_Along_Y_Dimension, Optical_Sectioning_Constant ) ;
                          
		if( Calibration_Along_X_Dimension > psf.NyqSpacialResolution() || 
		    Calibration_Along_Y_Dimension > psf.NyqSpacialResolution() )
		{
			double Radial_Calibration = psf.NyqSpacialResolution() ;
			if( Calibration_Along_X_Dimension < Radial_Calibration ) Radial_Calibration = Calibration_Along_X_Dimension ;
			if( Calibration_Along_Y_Dimension < Radial_Calibration ) Radial_Calibration = Calibration_Along_Y_Dimension ;
			Radial_Calibration /=  8.0 ;
			int    dimR = (int)( sqrt( Calibration_Along_X_Dimension * Calibration_Along_X_Dimension 
			                           * (dimx/2+1) * (dimx/2+1) + 
			                           Calibration_Along_Y_Dimension * Calibration_Along_Y_Dimension 
			                           * (dimy/2+1) * (dimy/2+1) ) 
			                    / Radial_Calibration ) * 2 ;                   		
			
			FluoRZPSF RZpsf( Objective_Numerical_Aperture, Emission_Light_Wavelength, 
			                 Act_Immersion_Medium_Refractive_Index, 
			                 dimR, dimz/2+1, Radial_Calibration, Optical_Sectioning_Constant ) ;
			if( mode == 1 || mode == 3 ) 
			{
				RZpsf.setImmersionMedium_Mismatch( Req_Immersion_Medium_Refractive_Index, 
				                                   Act_Immersion_Medium_Refractive_Index, 
				                                   Objective_Working_Distance ) ;
			}	
			if( mode == 2 || mode == 3 )
			{
				RZpsf.setCoverSlip_Mismatch( Req_Cover_Thickness, Act_Cover_Thickness, 
				                             Req_Cover_Refractive_Index, Act_Cover_Refractive_Index ) ;
			}			
			
			RZpsf.create( CheckPSFgeneration ) ;
			RZpsf.get3Dpsf( dimx, dimy, dimz, Calibration_Along_X_Dimension, Calibration_Along_Y_Dimension, 
			                Optical_Sectioning_Constant, MyPsf.data(), filename ) ;
		}
		else
		{
			if( mode == 1 || mode == 3 ) 
			{
				psf.setImmersionMedium_Mismatch( Req_Immersion_Medium_Refractive_Index, 
				                                 Act_Immersion_Medium_Refractive_Index, 
				                                 Objective_Working_Distance ) ;
			}			
			if( mode == 2 || mode == 3 )
			{
				psf.setCoverSlip_Mismatch( Req_Cover_Thickness, Act_Cover_Thickness, 
				                           Req_Cover_Refractive_Index, Act_Cover_Refractive_Index ) ;
			}
        		
        		psf.create( MyPsf.length(), MyPsf.width(), MyPsf.height(), MyPsf.data(), CheckPSFgeneration ) ;
			psf.exportProfile( filename ) ;
		}
		MyPsf.write( argv[1] ) ;
	}
	else
	{
		std::cout << usage ;
		return 1 ;
	} 	

	return 0 ;
}
