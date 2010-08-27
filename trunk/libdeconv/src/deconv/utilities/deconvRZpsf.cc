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
 * Filename:  deconvRZpsf.cc
 */


#include <stdlib.h>
#include "deconv/FluoRZPSF.h"


/*
 *	deconvRZpsf is designed for creating a 2-D RZpsf using the FluoRZPSF class.
 *
 *	Six input arguments are required by the program:
 *		<RZpsf> <PSF.para> <RidialDimension> <Sections> <RadialCalibration> <SectioningConstant>
 *
 *	The generated 2-D RZpsf table will be saved using FluoRZpsf::save():
 *		The header file : <RZpsf>.rzh
 *		The data   file : <RZpsf>.rzd 
 */


#define CheckPSFgeneration		true


char usage[] = "$deconvRZpsf <RZpsf> <PSF.para> <RadialDimension> <Sectons> <RadialCalibration> <SectioningConstant>\n" ;


int main( int argc, char ** argv )
{
	FluoRZPSF RZpsf ;
	char      line[1024] ;
	double    Objective_Numerical_Aperture, Emission_Light_Wavelength, Act_Immersion_Medium_Refractive_Index ; 
	double    Req_Immersion_Medium_Refractive_Index, Objective_Working_Distance ;
	double    Req_Cover_Refractive_Index, Act_Cover_Refractive_Index, Req_Cover_Thickness, Act_Cover_Thickness ;
	int       mode ;

	if( argc == 7 )
	{
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
		
		int    dimR = atoi( argv[3] ) ;
		int    dimZ = atoi( argv[4] ) ;
		double Radial_Calibration = atof( argv[5] ) ;
		double Optical_Sectioning_Constant = atof( argv[6] ) ;

		RZpsf.init( Objective_Numerical_Aperture, Emission_Light_Wavelength, Act_Immersion_Medium_Refractive_Index, 
		            dimR, dimZ, Radial_Calibration, Optical_Sectioning_Constant ) ;
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
		RZpsf.save( argv[1] ) ;
	}
	else
	{
		std::cout << usage ; 
		return 1 ;
	} 	

	return 0 ;
}
