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
 * Filename:  ViewerUtils.cc
 */


#include <math.h>
#include <stdio.h>
#include "deconv/MYdef.h"
#include "deconv/MYerror.h"
#include "ViewerUtils.h"


bool IsNotPowerOf2( int num )
{
	int newnum, pcnt ;
        
	pcnt = 0 ;
        newnum = num ;
        while( newnum > 0 ) 
        {
		newnum = newnum >> 1 ; 
		pcnt++ ; 
	}
        pcnt-- ;
	newnum = 1 << pcnt;

	if( newnum != num ) 
	{
		return true ;
	}
	else    return false ;
}


void name_seq_files( char * filehead, char * setname, int sliceID, int sliceEndID )
{
	char   sliceIDname[10], tempIDname[10] ;
	int    maxDigit, curDigit ;
	double logdigit ;

	if( sliceEndID < 10 )
	{
		sprintf( filehead, "%s%d", setname, sliceID) ;
	}
	else 
	{
		logdigit = (double) log10( (double) sliceEndID ) ;
		if( logdigit == ceil( logdigit ) ) maxDigit = 1 + (int) ceil( logdigit ) ;
		else                               maxDigit = (int) ceil( logdigit ) ;
               	
		logdigit = (double) log10( (double) sliceID ) ;
		if( logdigit == ceil( logdigit ) ) curDigit = 1 + (int) ceil( logdigit ) ;
		else                               curDigit = (int) ceil( logdigit ) ;
      		
		sprintf( sliceIDname, "%d", sliceID ) ;
		for( int i = curDigit ; i < maxDigit ; i++ ) 
		{
			sprintf( tempIDname, "%s", sliceIDname ) ;
			sprintf( sliceIDname, "0%s", tempIDname ) ;
		}
		sprintf( filehead, "%s%s", setname, sliceIDname );
	}
}


void write_my_byte_ppm( char* filehead, int width, int height, unsigned char * buf )
{
	std::string filename = std::string( filehead ) + ".ppm" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
	if( !fp ) 
	{
		throw ErrnoError( filename ) ;
	}
	else
	{
		fprintf( fp, "P6\n%d %d\n%d\n", width, height, 255 ) ;
		for( int i = height-1 ; i >=0 ; i-- )
		{
			if( fwrite( &buf[ i * width * 3 ], width, 3, fp ) != 3 ) throw WriteDataError( filename ) ;
		}
		fclose( fp ) ;
	}
}


unsigned char byte_contrast( double intensity, double IntensityWindow, double IntensityCenter )
{
	int temp = (int) floor( 0.5 + 255.0 / IntensityWindow * ( intensity - IntensityCenter + IntensityWindow/2.0 ) ) ;
	if      ( temp >= 255 )   return( ((unsigned char) 255)  ) ;
	else if ( temp <= 0   )   return( ((unsigned char) 0)    ) ;
	else                      return( ((unsigned char) temp) ) ;
}


unsigned short short_contrast( double intensity, double IntensityWindow, double IntensityCenter )
{
	int temp = (int) floor( 0.5 + 65535.0 / IntensityWindow * ( intensity - IntensityCenter + IntensityWindow/2.0 ) ) ;
	if      ( temp >= 65535 )   return ( ((unsigned short) 65535) ) ;
	else if ( temp <= 0     )   return ( ((unsigned short) 0)     ) ;
	else                        return ( ((unsigned short) temp)  ) ;
}


double contrast( double intensity, double IntensityWindow, double IntensityCenter )
{
	if( intensity >= (IntensityCenter+IntensityWindow/2.0) )      
	{
		return ( IntensityCenter+IntensityWindow/2.0 ) ;
	}
	else if( intensity <= (IntensityCenter-IntensityWindow/2.0) ) 
	{
		return ( IntensityCenter-IntensityWindow/2.0 ) ;
	}
	else 
	{
		return ( intensity ) ;
	}
}


void set_cylinder( int & rx, int & ry, int & hz )
{
	std::cout << " Input Cylinder Radius along X <<<< " ;
	std::cin  >> rx ;
	std::cout << " Input Cylinder Radius along Y <<<< " ;
	std::cin  >> ry ;
	std::cout << " Input Cylinder Height along Z <<<< " ;
	std::cin  >> hz ;
}


void set_ellipse( int & rx, int & ry, int & rz )
{
	std::cout << " Input Ellipse Radius along X <<<< " ;
	std::cin  >> rx ;
	std::cout << " Input Ellipse Radius along Y <<<< " ;
	std::cin  >> ry ;
	std::cout << " Input Ellipse Radius along Z <<<< " ;
	std::cin  >> rz ;
}
