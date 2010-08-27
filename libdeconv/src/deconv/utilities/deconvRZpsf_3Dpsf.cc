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
 * Filename:  deconvRZpsf_3Dpsf.cc
 */
 

#include <stdlib.h>
#include "deconv/FluoRZPSF.h"
#include "deconv/CCube.h"


/*
 *	deconvRZpsf_3Dpsf is designed for creating a 3-D PSF from a 2-D RZ PSF table using the FluoRZPSF and CCube classes.
 *
 *	Eight input arguments are required by the program:
 *		<RZpsf> <psf> <DimensionX> <DimensionY> <DimensionZ> <CalibrationX> <CalibrationY> <Sectioning Constant>
 *
 *	The 2-D RZpsf talbe profile will be read from <RZpsf>.rzh
 *	The 2-D RZpsf talbe data    will be read from <RZpsf>.rzd
 *
 *	The generated 3-D PSF will be saved using CCube::save():
 *		The header file : <psf>.hdr
 *		The data   file : <psf>.f32 if <DataType> is float. 
 *		                  <psf>.f64 if <DataType> is double.
 */


#define DataType           float                   
#define CheckPSFgeneration true


char usage[] = 
"$deconvRZpsf_3Dpsf <RZpsf> <psf> <DimensionX> < DimensionY> <DimensionZ> <CalibrationX> <CalibrationY> <Sectioning Constant>\n" ;


int main( int argc, char ** argv )
{
	FluoRZPSF        RZpsf ;
	CCube<DataType>  MyPsf ;
	char             filename[256] ;

	if( argc == 9 )
	{
		RZpsf.read( argv[1] ) ;
		sprintf( filename, "%s.log", argv[2] ) ;
		int    dimx = atoi( argv[3] ) ;
		int    dimy = atoi( argv[4] ) ;
		int    dimz = atoi( argv[5] ) ;
		double Calibration_Along_X_Dimension = atof( argv[6] ) ; 
		double Calibration_Along_Y_Dimension = atof( argv[7] ) ;
		double Optical_Sectioning_Constant   = atof( argv[8] ) ;
		
		MyPsf.init( dimx, dimy, dimz ) ;

		if( RZpsf.get3Dpsf( dimx, dimy, dimz, Calibration_Along_X_Dimension, Calibration_Along_Y_Dimension, 
		                    Optical_Sectioning_Constant, MyPsf.data(), filename ) )
		{ 
			std::cout << " program failed to create the 3-D PSF.\n" ; 
			return 1 ;
		}
		else
		{
			MyPsf.write( argv[2] ) ;
			return 0 ;
		}
	}
	else
	{
		std::cout << usage ; 
		return 1 ;
	}
}
