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
 * Filename:  FluoRZPSF.h
 */
 

#ifndef FLUORZPSF_H
#define FLUORZPSF_H


#include "FluoPSF.h"


/*
 *	======================================================================================================
 *	FluoRZPSF is used to create a 2-D RZ PSF given physical microscope parameters and set up based on 
 *	a theoretical PSF Model developed by Gibson and Lanni at CMU. The 2-D RZ PSF will be used to 
 *	simulate a 3-D PSF volume for deconvolution. Please see the description for a 3-D PSF in "Fluo3DPSF.h"
 *	======================================================================================================
 *
 *	DimR (R) is the slow varying dimension of a RZ PSF and indexes from r = 0 to RadialDimension()-1.
 *	It is equivalent to the DimensionX and DimensionY of a 3-D PSF due to the 3-D PSF is cylindrically symmetric.
 *
 *	Sections (S) is the fast varying dimension of a RZ PSF and indexes from s = 0 to OpticalSecitons()-1.
 *	It is equivalent to the DimensionZ of a 3-D PSF.
 *	
 *	The RZ PSF data is stored in the FluoRZPSF one-dimensional <DoubleArray> as "r + s * RadialDimension()".
 *
 *	A RZ PSF is saved in two files: a header file and a data file.
 *	The header file is a text file and its suffix is ".rzh".
 *	The data file is a binary file storing the RZpsf data (one-dimensional array) in a binary stream. 
 *
 *
 *	=============================================================================================================
 *	How to generate a 3-D PSF : 3Dpsf(X,Y,Z) from a 2-D RZ PSF : RZpsf(R,S) under the following conditions:
 *	3Dpsf calibrations = (dx, dy, dz); RZpsf calibrations = (dr, dzz); the compensation parameter Points2Sum = P.
 *	=============================================================================================================
 *
 *		------------------------------
 *		Requirement for the parameters
 *		------------------------------
 *			P is an odd positive integer.
 *			dz = dzz * N where N is an positive integer, and dx >= dr * P, and dy >= dr * P. 
 *			X < R * dr / dx / 2, and  Y < R * dr / dy / 2, and Z < S * dzz / dz.
 *
 *		-----------------------------------------------------
 *		Calculation of an individual 3Dpsf voxel at (x, y, z)
 *		----------------------------------------------------- 
 *			It is equivalent to calculate RZpsf(r = sqrt(x^2+y^2), z) by performing RZ PSF lookups
 *			and spline interpolation on an one-dimensional array RZpsf(0:R-1, z*dz/dzz) in the table.
 *
 *		----------------------------------------------------
 *		Application of the compensation using Points2Sum = P
 *		---------------------------------------------------- 
 *			For each 3Dpsf plane, it is first to calculate a plane whose size is X*P by Y*P, and then
 *			to generate the X by Y plane of the 3Dpsf by summing P*P pixels on the X*P by Y*P plane. 
 */


class RZHError1 : public Error
{
	public: 
	RZHError1( char * filename )
	{
		_error << " It is a bad .rzh file which does not have enough information : " << filename << "\n" ;		
	}
} ;


class RZHError2 : public Error
{
	public: 
	RZHError2( char * filename )
	{
		_error << " It is a bad .rzh file whose Nyquist samplings do not match with its NA and WL : " 
		       << filename << "\n" ;
	}
} ;


class FluoRZPSF : public FluoPSF 
{
	public:
	virtual ~FluoRZPSF() {}
	
	
	/*
	 *	Default Constructor
	 */
	FluoRZPSF() { _NA = 0.0 ; }
	
	
	/*
	 *	Constructor
	 *	Input:
	 *		NA,                 it is the objective numberical aperture.
	 *		WL,                 it is the emission light wavelength.
	 *		RI,                 it is the refractive index of the immersion medium used for the objective.
	 *		RadialCalibration,  it is the radial sampling rate (in microns) in the object space.
	 *		SectioningConstant, it is the sectioning rate along the optical axis (in microns).
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 */
	FluoRZPSF( double NA, double WL, double RI, int DimR, int Sections, double RadialCalibration, double SectioningConstant ) ;
	
	
	/*
	 *	Constructor
	 *	Input:
	 *		NA,                 it is the objective numberical aperture.
	 *		WL,                 it is the emission light wavelength.
	 *		RI,                 it is the refractive index of immersion medium used for the objective.
	 *		DimR,               it is the fast varying dimension of a RZ PSF.
	 *		DimZ,               it is the slow varying dimension of a RZ PSF.
	 *		RadialCalibration,  it is the radial sampling rate (in microns) of a RZ PSF.
	 *		SectioningConstant, it is the sectioning rate (in microns).
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 */
	void init( double NA, double WL, double RI, int DimR, int Sections, double RadialCalibration, double SectioningConstant ) ;
	
	
	/*
	 *	Get pravite members
	 *	RadialCalibration() returns the radial sampling rate (in microns) of a RZ PSF.
	 *	RadialDimension()   returns the fast varying dimension of a RZ PSF.
	 *	OpticalSecions()    returns the slow varying dimension of a RZ PSF.
	 *	RZpsfData()         returns the double array storing the data of a RZ PSF.
	 */
	double 	RadialCalibration()  { return _DR ;          }
	int    	RadialDimension()    { return _DimR ;        }
	int     OpticalSecions()     { return _Sections;     }
	double* RZpsfData()          { return _RZpsf ; }


	/*
	 *	Calculate maximum dimensions of a 3-D PSF that a RZ PSF can generate given its calibrations 
	 *	Input:
	 *		CalibraionX,        it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel width (the fast varying dimension of the image).
	 *		CalibraionY,        it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel height (the slow varying dimension of the image).
	 *		SectioningConstant, it is the distance (in microns) between two successive sectioning images.
	 *	Output:
	 *		return the maximum value of a dimension. 
	 */
	int     maxDimensionX  ( double CalibrationX ) ;
	int     maxDimensionY  ( double CalibrationY ) ;
	int     maxSections    ( double SectioningConstant ) ;
	
	
	/*
	 *	Calculate minimum lateral calibration of a 3-D PSF that a RZ PSF can generate given a compensation parameter 
	 *	Input:  
	 *		Points2Sum, it is the compensation parameter.
	 *	Output: 
	 *		return the minimum lateral calibration.
	 */
	double 	minCalibration( int Points2Sum ) ;
	
	
	/*
	 *	Save RZ PSF data to a pair of header (.rzh) and data (.rzd) files
	 *	Input:
	 *		filehead, it is the file name without suffix.
	 *	Throw:
	 *		throw an error if fail.
	 */
	void 	save( const char * filehead ) ;
	
	
	
	/*
	 *	Read RZ PSF data from a pair of header (.rzh) and data (.rzd) files
	 *	Input:
	 *		filehead, it is the file name without suffix.
	 *	Throw:
	 *		throw an error if fail.
	 */
	void 	read( const char * filehead ) ;
	
	
	/*
	 *	Create a 2-D RZ PSF
	 *	Input:
	 *		Check, it is to indicate whether to print the generation progress or not; its default is true. 
	 *	Throw:
	 *		throw an error if the RZ PSF has not been initialized.
	 */
	void 	create( bool Check = true ) ;
	
	
	/*
	 *	Generate a 3-D PSF in the single/double floating precision from a 2-D RZ PSF
	 *	Input:
	 *		nx,         it is the DimensionX of a 3-D PSF (the fastest varying dimension).
	 *		ny,         it is the DimensionY of a 3-D PSF.
	 *		nz,         it is the DimensionZ of a 3-D PSF (the slowest varying dimension).
	 *		dx,         it is the corresponding length (in microns) in the object space to
	 *		            a sectioning image pixel width (the fast varying dimension of the image).
	 *		dy,         it is the corresponding length (in microns) in the object space to
	 *		            a sectioning image pixel height (the slow varying dimension of the image).
	 *		dz,         it is the distance (in microns) between two successive sectioning images.
	 *		psf,        it is an one-dimensional nx*ny*nz array holding the 3-D PSF data.
	 *		file,       it is the name of the exported profile for the generated 3-D PSF; 
	 *		            its default value is NULL in which case the profile will not be exported. 
	 *		Points2Sum, it is the compensation parameter; its default value is 5.
	 *	Output:
	 *		return 0 if success.
	 *		return 1 if given an invalid DimensionX (nx).
	 *		return 2 if given an invalid DimensionY (ny).
	 *		return 3 if given an invalid DimensionZ (nz).
	 *		return 4 if given an invalid compensation parameter (Points2Sum).
	 *		return 5 if given an invalid X-calibration (dx) or Y-calibration (dy).
	 *		return 6 if given an invalid sectioning constant (dz).	
	 *	Throw:
	 *		throw an error if the RZ PSF has not been initialized.	
	 */	 
	int  	get3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, float  * psf, const char * file = NULL, int Points2Sum = 5 ) ;
	int  	get3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, double * psf, const char * file = NULL, int Points2Sum = 5 ) ;
	
		
	private:
	int         _DimR, _Sections ;
	double      _DR ;
	double*     _RZpsf ;
	int         _check3Dpsf( int nx, int ny, int nz, double dx, double dy, double dz, int Points2Sum ) ;
} ;


#endif    /*   #include "FluoRZPSF.h"   */
