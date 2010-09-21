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
 * Filename:  FluoPSF.cc
 */
 
#include <stdio.h>
#include <math.h>
#include "FluoPSF.h" 


double FluoPSF_OPD( double rho, FluoPSF_func_params * p )
{
	double temp = sqrt( p->ActImmRI2 - rho ) ;
	double  opd = p->defocus * ( temp - p->ActImmRI ) ;
		
	if( fabs( p->ActImmRI - p->ReqImmRI ) > DiffEpsilon )
	{
		opd += ( p->WD * ( sqrt(p->ReqImmRI2 - rho) - p->ActImmRI / p->ReqImmRI * temp - 
		                                    p->ReqImmRI + p->ActImmRI2 / p->ReqImmRI ) ) ;
	}
		
	if( ( fabs( p->ActCovTh - p->ReqCovTh ) > DiffEpsilon || fabs( p->ActCovRI - p->ReqCovRI ) > DiffEpsilon ) 
	      && p->ReqCovTh2 > rho && p->ActCovTh2 > rho )    
	{
		opd += ( p->ReqCovTh * ( sqrt(p->ReqCovTh2 - rho) - p->ActImmRI / p->ReqCovRI * temp ) -
		         p->ActCovTh * ( sqrt(p->ActCovTh2 - rho) - p->ActImmRI / p->ActCovRI * temp ) - 
		         p->CovThxRI + p->ActImmRI2 * p->CovTh$RI ) ;
	}
	
	return opd ;
}



double FluoPSF_cos_func( double rho, void * params )
{
	struct FluoPSF_func_params * p = (struct FluoPSF_func_params *) params ;
	
	if( p->ActImmRI2 > rho ) 
	{
		double opd = FluoPSF_OPD( rho, p ) ;
		return ( j0( p->WN * p->radius * sqrt(rho) ) * cos( p->WN * opd ) ) ;
	}
	else
	{
		return 0.0 ;
	}
}



double FluoPSF_sin_func( double rho, void * params )
{
	struct FluoPSF_func_params * p = (struct FluoPSF_func_params *) params ;
	  
	if( p->ActImmRI2 > rho ) 
	{
		double opd = FluoPSF_OPD( rho, p ) ;
		return ( j0( p->WN * p->radius * sqrt(rho) ) * sin( p->WN * opd ) ) ;
	}
	else
	{
		return 0.0 ;
	}
}



/* public functions */

void FluoPSF::setImmersionMedium_Mismatch( double RequiredRI, double ActualRI, double WD )
{
	if( _NA > 0.0 )
	{
		if( RequiredRI >= RI_LowerLimit && RequiredRI <= RI_UpperLimit )
		{
			 _ReqImmRI = RequiredRI ;
			__ReqImmRI = _ReqImmRI * _ReqImmRI ;
		}
		else	throw RIError( RequiredRI ) ;
		
		if( ActualRI >= RI_LowerLimit && ActualRI <= RI_UpperLimit ) 
		{
			 _ActImmRI = ActualRI ;
			__ActImmRI = _ActImmRI * _ActImmRI ;
		}
		else	throw RIError( ActualRI ) ;
	
		if( WD >= WD_LowerLimit && WD <= WD_UpperLimit )
		{
			_WD = WD ;
		}
		else	throw WDError( WD ) ;
	}
	else 	throw PSFError( 0 ) ;	
}



void FluoPSF::setCoverSlip_Mismatch( double RequiredTh, double ActualTh, double RequiredRI, double ActualRI )
{
	if( _NA > 0.0 )
	{
		if( RequiredTh >= Cover_LowerLimit && RequiredTh <= Cover_UpperLimit )
		{
			 _ReqCovTh = RequiredTh ;
			__ReqCovTh = _ReqCovTh * _ReqCovTh ;
		}
		else	throw CoverError( RequiredTh ) ;
		
		if( ActualTh >= Cover_LowerLimit && ActualTh <= Cover_UpperLimit )
		{
		 	_ActCovTh = ActualTh ;
			__ActCovTh = _ActCovTh * _ActCovTh ;
		}
		else 	throw CoverError( ActualTh ) ;
		
		if( RequiredRI >= RI_LowerLimit && RequiredRI <= RI_UpperLimit )
		{
			 _ReqCovRI = RequiredRI ;
		}
		else	throw RIError( RequiredRI ) ;
	
		if( ActualRI >= RI_LowerLimit &&   ActualRI <= RI_UpperLimit ) 
		{
		 	_ActCovRI = ActualRI ;
		}
		else	throw RIError( ActualRI ) ;
	
		_CovThxRI = _ReqCovTh * _ReqCovRI - _ActCovTh * _ActCovRI ;
		_CovTh$RI = _ReqCovTh / _ReqCovRI - _ActCovTh / _ActCovRI ;
	}
	else	throw PSFError( 0 ) ;
}
 


/* protected functions */

void FluoPSF::_init( double NA, double WL, double RI )
{
	if( NA >= NA_LowerLimit && NA <= NA_UpperLimit ) 
	{
		_NA = NA ;
	}
	else	throw NAError( NA ) ;
	
	if( WL >= WL_LowerLimit && WL <= WL_UpperLimit ) 	
	{
		_WL = WL ;
		_WN = 2.0 * M_PI / _WL ;
	}
	else	throw WLError( WL ) ;
	
	if( RI >= RI_LowerLimit && RI <= RI_UpperLimit )
	{
		_ActImmRI  = _ReqImmRI  = RI ;
		__ActImmRI = __ReqImmRI = (RI * RI) ;
	}
	else	throw RIError( RI ) ;
	
	_NyqDXY    = WL / ( 4.0 * NA ) ;
	_NyqDF     = WL / (  NA * NA ) ;
	_DZ        = 0.0 ;
	_WD        = 0.0 ;
	_ReqCovTh  = 0.0 ;
	__ReqCovTh = 0.0 ;
	_ActCovTh  = 0.0 ;
	__ActCovTh = 0.0 ;
	_ReqCovRI  = 0.0 ;
	_ActCovRI  = 0.0 ;
	_CovThxRI  = 0.0 ;
	_CovTh$RI  = 0.0 ;
}



void FluoPSF::_exportCommon( FILE * fp )
{
	fprintf( fp, "%10.4f -> Nyquist Spacial   Resolution (um) at the focus plane in the object space.\n", _NyqDXY ) ;
	fprintf( fp, "%10.4f -> Nyquist Depth of Focus Field (um) along the optical axis in the object space.\n", _NyqDF ) ;
	fprintf( fp, "\n" ) ;
	
	fprintf( fp, "%10.4f -> Microscope Objective Numerical Aperture.\n", _NA ) ;
	fprintf( fp, "%10.4f -> Fluorescence Emission Light Wavelength (um).\n", _WL ) ;
	fprintf( fp, "%10.4f -> Refractive Index of Immersion Medium used for Microscope Objective.\n", _ActImmRI ) ;
	fprintf( fp, "\n" ) ;
	
	if( fabs( _ActCovRI - _ReqCovRI ) > DiffEpsilon || fabs( _ActCovTh - _ReqCovTh ) > DiffEpsilon ) 
	{
		fprintf( fp, "%10.4f -> Refractive Index of Cover Slip/Glass used to contain specimens.\n", _ActCovRI ) ;
		fprintf( fp, "%10.4f -> Refractive Index of Cover Slip/Glass required by Microscope Objective.\n", _ReqCovRI ) ;
		fprintf( fp, "%10.4f -> Thickness   (um) of Cover Slip/Glass used to contain specimens.\n", _ActCovTh ) ;
		fprintf( fp, "%10.4f -> Thickness   (um) of Cover Slip/Glass required by Microscope Objective.\n", _ReqCovTh ) ;
		fprintf( fp, "\n" ) ;
	}

	if( fabs( _ActImmRI - _ReqImmRI ) > DiffEpsilon ) 
	{
		fprintf( fp, "%10.4f -> Refractive Index of Immersion Medium required by Microscope Objective.\n", _ActImmRI ) ;
		fprintf( fp, "%10.4f -> Microscope Objective Working Distance (um).\n", _WD ) ;
	}	
	fprintf( fp, "\n" ) ;
}



double FluoPSF::_IntegralPSF( gsl_integration_workspace * ws, double lolimit, double uplimit, double r, double defocus )
{
	gsl_function                cos_func, sin_func ;
	struct FluoPSF_func_params  params = { r, defocus, _WN, _WD, _ActImmRI, __ActImmRI, _ReqImmRI, __ReqImmRI, _ReqCovTh, 
	                                       __ReqCovTh, _ActCovTh, __ActCovTh, _ReqCovRI, _ActCovRI, _CovThxRI, _CovTh$RI  } ;
	double                      cosret, sinret, coserr, sinerr ;
		   
	cos_func.function = &FluoPSF_cos_func ;
	cos_func.params   = &params ;
	sin_func.function = &FluoPSF_sin_func ;
	sin_func.params   = &params ;
	
	gsl_integration_qag( &cos_func, lolimit, uplimit, IntegrationEpsabs, IntegrationEpsrel, IntegrationLimit, IntegrationKey, ws, &cosret, &coserr ) ;
	gsl_integration_qag( &sin_func, lolimit, uplimit, IntegrationEpsabs, IntegrationEpsrel, IntegrationLimit, IntegrationKey, ws, &sinret, &sinerr ) ;
			      
	return ( cosret * cosret + sinret * sinret ) ;
}
