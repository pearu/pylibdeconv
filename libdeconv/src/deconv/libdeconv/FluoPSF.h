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
 * Filename:  FluoPSF.h
 */
 
 
#ifndef FLUOPSF_H
#define FLUOPSF_H


#include <gsl/gsl_integration.h>
#include "MYdef.h"
#include "MYerror.h"


#define DiffEpsilon            1.0E-10  // two physical parameters are regarded to be same if
                                        // their absolute difference is more than DiffEpsilon.  

#define IntegrationEpsabs      1.0E-9   // desired absolute error limit
#define IntegrationEpsrel      1.0E-6   // desired relative error limit
#define IntegrationKey         6        // 61-point Gauss-Kronrod rule 
#define IntegrationLimit       100      // maximum number of subintervals 

#define NA_LowerLimit          0.2
#define NA_UpperLimit          2.0

#define WL_LowerLimit          0.1      // um
#define WL_UpperLimit          1.0      // um

#define WD_LowerLimit          10.0     // um
#define WD_UpperLimit          15000.0  // um

#define RI_LowerLimit 	       0.5
#define RI_UpperLimit          2.0

#define Cover_LowerLimit       10.0     // um
#define Cover_UpperLimit       10000.0  // um

#define Calibration_LowerLimit 0.001    // um
#define Calibration_UpperLimit 10.00    // um

#define Sectioning_LowerLimit  0.01     // um
#define Sectioning_UpperLimit  100.0    // um

#define Samples_LowerLimit     8


/*
 *	Structure of physical pamameters used to calculate an integration point
 */
struct FluoPSF_func_params 
{ 
	double radius ;     // radius
	double defocus ;    // defous
	double WN ;         // wave number
	double WD ;         // working distance
	double ActImmRI ;   // refractive index of the immersion medium used for the objective
	double ActImmRI2 ;  // ActImmRI * ActImmRI
	double ReqImmRI ;   // refractive index of the immersion medium required by the objective
	double ReqImmRI2 ;  // ReqImmRI * ReqImmRI
	double ReqCovTh ;   // thickness of the cover slip/glass required by the objective
	double ReqCovTh2 ;  // ReqCovTh * ReqCovTh
	double ActCovTh ;   // thickness of the cover slip/glass used for the objective
	double ActCovTh2 ;  // ActCovTh * ActCovTh
	double ReqCovRI ;   // refractive index of the cover slip/glass required by the objective
	double ActCovRI ;   // refractive index of the cover slip/glass used for the objective
	double CovThxRI ;   // _ReqCovTh * _ReqCovRI - _ActCovTh * _ActCovRI 
	double CovTh$RI ;   // _ReqCovTh / _ReqCovRI - _ActCovTh / _ActCovRI
} ;


/*
 *	Calcuate the optical path difference (OPD) given an integration point
 *	Input:
 *		rho, it is the integral variable. 
 *		p,   it is the structure of physical parameters used to calculate an integration point.
 *	Output:
 *		the OPD value at the given integration point.
 */
double FluoPSF_OPD( double rho, FluoPSF_func_params * p ) ;


/*
 *	Calcuate the cosine value at an integration point
 *	Input:
 *		rho, it is the integral variable. 
 *		p,   it is the structure of physical parameters used to calculate an integration point.
 *	Output:
 *		the cosine value at the given integration point.
 */
double FluoPSF_cos_func( double rho, void * params ) ; 


/*
 *	Calcuate the sine value at an integration point
 *	Input:
 *		rho, it is the integral variable. 
 *		p,   it is the structure of physical parameters used to calculate an integration point.
 *	Output:
 *		the sine value at the given integration point.
 */
double FluoPSF_sin_func( double rho, void * params ) ;


class NAError : public Error
{
	public:
	NAError( double NA )
	{
		_error << " Objective Numerical Aperture = " << NA << " must be set within " 
		       << NA_LowerLimit << " ~ " << NA_UpperLimit << ".\n" ;
	}
} ;


class WLError : public Error
{
	public:
	WLError( double WL )
	{
		_error << " Emission Light Wavelength = " << WL << " (um) must be set within " 
		       << WL_LowerLimit << " ~ " << WL_UpperLimit << ".\n" ;
	}
} ;


class WDError : public Error
{
	public:
	WDError( double WD )
	{
		_error << " Objective Working Distance = " << WD << " (um) must be set within " 
		       << WD_LowerLimit << " ~ " << WD_UpperLimit << ".\n" ;
	}
} ;


class RIError : public Error
{
	public:
	RIError( double RI )
	{
		_error << " Refractive Index = " << RI << " must be set within " 
		       << RI_LowerLimit << " ~ " << RI_UpperLimit << ".\n" ;
	}
} ;


class CoverError : public Error
{
	public:
	CoverError( double Th )
	{
		_error << " Cover Slip/Glass Thickness = " << Th << " (um) must be set within " 
		       << Cover_LowerLimit << " ~ " << Cover_UpperLimit << ".\n" ;
	}
} ;


class CalibrationError : public Error
{
	public:
	CalibrationError( double c )
	{
		_error << " X/Y Calibration = " << c << " (um/pixel) in the object space must be set within " 
		       << Calibration_LowerLimit << " ~ " << Calibration_UpperLimit << ".\n" ;
	}
} ;


class SectioningError : public Error
{
	public:
	SectioningError( double DZ )
	{
		_error << " Sectioning Constant = " << DZ << " (um/step) must be set within " 
		       << Sectioning_LowerLimit << " ~ " << Sectioning_UpperLimit << ".\n" ;
	}
} ;


class SamplesError : public Error
{
	public:
	SamplesError( int sample )
	{
		_error << " DimensionX/Y/R/Z Samples = " << sample << " mustNOT be smaller than "<< Samples_LowerLimit << ".\n" ;
	}
} ;


class PSFError : public Error
{
	public:
	PSFError( int type )
	{
		switch( type )
		{
			case 0:
				_error << " Fluo3DPSF/FluoRZPSF hasNOT been initialized yet : use init().\n" ;
				break ;
				
			case 1:
				_error << " Fluo3DPSF hasNOT been initialized yet : use init().\n" ;
				break ;
				
			case 2:
				_error << " FluoRZPSF hasNOT been initialized yet : use init().\n" ;
				break ;
				
			case 3:
				_error << " FluoRZPSF table data hasNOT been created yet : use create().\n" ;
				break ;
				
			default:
				break ;
		}
	}
};


/*
 *	The FluoPSF is a base class for the Fluo3DPSF and FluoRZPSF classes.
 */


class FluoPSF 
{
	public:
	virtual ~FluoPSF() {}
	
	
	/*
	 *	Get protected members
	 *	NumbericalAperture()   returns the objective numerical aperture (NA).
	 *	EmissionWavelength()   returns the emission light wavelength (WL) in microns.
	 *	EmissionWaveNumber()   returns the emission light wave-number and equals to WL/2/pi.
	 *	ImmersionMediumRI()    returns the refractive index of immersion medium used for the objective.
	 *	NyqSpacialMediumRI()   returns the Nyquist sampling rate (in microns) calculated from NA and WL. 
	 *	NyqDepthOfFocusField() returns the Nyquist depth of field (in microns) calculated from NA and WL.
	 *	SectioningConstant()   returns the optical sectioning constant (in microns) along the optical axis. 
	 */
	double NumericalAperture()    { return _NA ;       }
	double EmissionWaveLength()   { return _WL ;       }
	double EmissionWaveNumber()   { return _WN ;       }
	double ImmersionMediumRI()    { return _ActImmRI ; }
	double NyqSpacialResolution() { return _NyqDXY ;   }
	double NyqDepthOfFocusField() { return _NyqDF ;    }
	double SectioningConstant()   { return _DZ ;       }
	
	
	/*
	 *	Set the improper use of the immersion medium for the objective
	 *	Input:
	 *		RequiredRI, it is the refractive index of the immersion medium required by the objective.
	 *		ActualRI,   it is refractive index of the immersion medium used for the objective.
	 *		WD,         it is the objective working distance (in microns).
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 *		throw an error if the objective numerical aperture is not larger than 0.
	 */
	void setImmersionMedium_Mismatch( double RequiredRI, double ActualRI, double WD ) ;
	
	
	/*
	 *	Get the imporper use of the immersion medium for the objective
	 *	Input:
	 *		RequiredRI, it is to store the refractive index of immersion medium required by the objective.
	 *		ActualRI,   it is to store the refractive index of immersion medium used for the objective.
	 *		WD,         it is to store the objective working distance (in microns).
	 */
	void getImmersionMedium_MisMatch( double & RequiredRI, double & ActualRI, double & WD )
	{
		RequiredRI = _ReqImmRI ;
		ActualRI   = _ActImmRI ;
		WD         = _WD ; 
	}
	
	
	/*
	 *	Set the improper use of the cover slip/glass for the objective
	 *	Input:
	 *		RequiredTh, it is the thickness (in microns) of the cover slip/glass required by the objective.
	 *		ActualTh,   it is the thickness (in microns) of the cover slip/glass used for the objective.
	 *		RequiredRI, it is the refractive index of the cover slip/glass required by the objective.
	 *		ActualRI,   it is the refractive index of the cover slip/glass used for the objective.
	 *	Throw:
	 *		throw an error if an input is out of the pre-defined range.
	 *		throw an error if the objective numerical aperture is not larger than 0.
	 */
	void setCoverSlip_Mismatch( double RequiredTh, double ActualTh, double RequiredRI, double ActualRI ) ;
	
	
	/*
	 *	Get the improper use of the cover slip/glass for the objective
	 *	Input:
	 *		RequiredTh, it is to store the thickness (in microns) of the cover slip/glass required by the objective.
	 *		ActualTh,   it is to store the thickness (in microns) of the cover slip/glass used for the objective.
	 *		RequiredRI, it is to store the refractive index of the cover slip/glass required by the objective.
	 *		ActualRI,   it is to store the refractive index of the cover slip/glass used for the objective.
	 */
	void getCoverSlip_Mismatch( double & RequiredTh, double & ActualTh, double & RequiredRI, double & ActualRI )
	{
		RequiredTh = _ReqCovTh ;
		ActualTh   = _ActCovTh ;
		RequiredRI = _ReqCovRI ;
		ActualRI   = _ActCovRI ;
	}
	
	
	protected:
	double _WN, _NA, _WL, _WD, _NyqDXY, _NyqDF , _DZ ;
	double _ActImmRI, __ActImmRI, _ReqImmRI, __ReqImmRI ;
	double _ReqCovTh, __ReqCovTh, _ActCovTh, __ActCovTh, _ReqCovRI, _ActCovRI, _CovThxRI, _CovTh$RI ;
	
	void   _init( double NA, double WL, double RI ) ;
	void   _exportCommon( FILE * fp ) ;
	double _IntegralPSF( gsl_integration_workspace * ws, double lolimit, double uplimit, double r, double defocus ) ;
} ;


#endif    /*   #include "FluoPSF.h"   */
