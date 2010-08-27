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
 * Filename:  Fluo3DPSF.h
 */
 

#ifndef FLUO3DPSF_H
#define FLUO3DPSF_H


#include "FluoPSF.h"


/*
 *	====================================================================================
 *	Fluo3DPSF is used to create a 3-D PSF volume given physical microscope parameters
 *	and set up based on a theoretical PSF Model developed by Gibson and Lanni at CMU.
 *	The simulated PSF will be used for deconvolution applied to fluorescence microscopy.
 *	====================================================================================
 *
 *	DimensionX (X) is the fastest varying dimension of a 3-D PSF and indexes from x = 0 to DimensionX()-1.
 *	It is corresponding to the fast varying dimension of a sectioning image.
 *
 *	DimensionY (Y) is the middle dimension of a 3-D PSF and indexes from y = 0 to DimensionY()-1.
 *	It is corresponding to the slow varying dimension of a sectioning image.
 *
 *	DimensionZ (Z) is the slowest varying dimension of a 3-D PSF and indexes from z = 0 to DimensionZ()-1.
 *	It is corresponding to the sectioning dimension along the optical axis.
 *	
 *	The generated 3-D psf data is stored in an one-dimensional array as "x+y*DimensionX()+z*DimensionX()*DimensionY()".
 */


class Fluo3DPSF : public FluoPSF 
{
	public:
	virtual ~Fluo3DPSF() {}
	
	/*
	 *	Default Constructor
	 */
	Fluo3DPSF() { _NA = 0.0 ; }
	
	
	/*
	 *	Constructor
	 *	Input:
	 *		NA,                 it is the objective numberical aperture.
	 *		WL,                 it is the emission light wavelength.
	 *		RI,                 it is the refractive index of immersion medium used for the objective.
	 *		CalibrationX,       it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel width (the fast varying dimension of the image).
	 *		CalibrationY,       it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel height (the slow varying dimension of the image).
	 *		SectioningConstant, it is the distance (in microns) between two successive sectioning images.
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 */
	Fluo3DPSF( double NA, double WL, double RI, double CalibrationX, double CalibrationY, double SectioningConstant ) ;


	/*
	 *	Set up physical parameters for a 3-D PSF 
	 *	Input:
	 *		NA,                 it is the objective numberical aperture.
	 *		WL,                 it is the emission light wavelength.
	 *		RI,                 it is the refractive index of immersion medium used for the objective.
	 *		CalibrationX,       it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel width (the fast varying dimension of the image).
	 *		CalibrationY,       it is the corresponding length (in microns) in the object space to
	 *		                    a sectioning image pixel height (the slow varying dimension of the image).
	 *		SectioningConstant, it is the distance (in microns) between two successive sectioning images.
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 */	
	void 	init( double NA, double WL, double RI, double CalibrationX, double CalibrationY, double SectioningConstant ) ;
	
	
	/*
	 *	Get pravite members
	 *	CalibrationX() returns the corresponding length (in microns) in the object space to
	 *	               a sectioning image pixel width (the fast varying dimension of the image).
	 *	CalibrationY() returns the corresponding length (in microns) in the object space to
	 *	               a sectioning image pixel height (the slow varying dimension of the image).
	 *	DimensionX()   returns the fastest varying dimension of a 3-D PSF.
	 *	DimensionY()   returns the middle          dimension of a 3-D PSF.
	 *	DimensionZ()   returns the slowest varying dimension of a 3-D PSF. 
	 */
	double 	CalibrationX()  { return _DX ;    }
	double 	CalibrationY()  { return _DY ;    }
	int    	DimensionX()    { return _DimX ;  }
	int    	DimensionY()    { return _DimY ;  }
	int     DimensionZ()    { return _DimZ ;  }
	
	
	/*
	 *	Create a 3-D PSF in the single/double precision given its dimensions :
	 *	Input:
	 *		DimX,  it is the fastest varying dimension of a 3-D PSF.
	 *		DimY,  it is the middle          dimension of a 3-D PSF.
	 *		DimZ,  it is the slowest varying dimension of a 3-D PSF. 
	 *		psf,   it is an one-dimensional DimX*DimY*DimZ array holding the PSF data.
	 *		Check, it is to indicate whether to print the generation progress or not; its default is true.
	 *	Throw:
	 *		throw an error if a dimension is not even or out of the pre-defined range.
	 *		throw an error if the objective numerical aperture is not larger than 0. 
	 */
	void 	create( int DimX, int DimY, int DimZ, float  * psf, bool Check = true ) ;
	void 	create( int DimX, int DimY, int DimZ, double * psf, bool Check = true ) ;
	
	
	/*
	 *	Export the profile of a 3-D PSF to a text file
	 *	Input:
	 *		filename, it is the name of the text file to be written including suffix. 
	 *	Throw:
	 *		throw an error if fail.
	 */
	void 	exportProfile( const char * filename ) ;
	
	private:
	int    _DimX, _DimY, _DimZ ;
	double _DX, _DY ;
	
	void   _setDimensions( int DimX, int DimY, int DimZ ) ; 
} ;


#endif    /*   #include "Fluo3DPSF.h"   */
