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
 * Filename:  MYcube.cc
 */


#include <stdio.h>
#include <unistd.h>
#include "MYcube.h"

void read_my_cube_hdr( std::string filehead, int& length, int& width, int& height )
{
	std::string filename = filehead + ".hdr" ;
	FILE * fp = fopen( filename.c_str(), "rt" ) ;
  	
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else
	{ 
		int h, w, l ;
		
		if (fscanf (fp, "%d %d %d", &h, &w, &l) != 3) 
			throw ReadDataError( std::string(filename) ) ;

		width = w ;
		height = h ;
		length = l ;
		
		fclose (fp) ;
	}
}

void read_my_byte_cube_data( const char* filename, int cube_size, unsigned char * buf )
{
	FILE * fp = fopen( filename, "r" ) ;
	
	if( !fp ) 
		throw ErrnoError( std::string(filename) ) ;
	else
	{
		int index = 0 ;
		unsigned char* pbuf = buf ;

		while (index < cube_size) {
			if (!fread (pbuf++, sizeof (unsigned char), 1, fp))
				throw ReadDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	} 
}

void read_my_short_cube_data( const char* filename, int cube_size, unsigned short * buf )
{
	FILE * fp = fopen( filename, "r" ) ;
	
	if( !fp ) 
		throw ErrnoError( std::string(filename) ) ;
	else
	{
		int index = 0 ;
		unsigned short* pbuf = buf ;

		while (!feof(fp)) {
			unsigned short val ;
			if (!fread (&val, sizeof (unsigned short), 1, fp))
				throw ReadDataError (std::string(filename)) ;
				
         	swab ((char*)&val, (char*)pbuf++, sizeof (unsigned short)) ;
			index ++ ;
		}

		fclose( fp ) ;

		if (index != cube_size)
			throw ReadDataError (std::string(filename)) ;		
	}
}

void read_my_single_cube_data( const char* filename, int cube_size, float * buf )
{
	FILE * fp = fopen( filename, "r" ) ;
	
	if( !fp ) 
		throw ErrnoError( std::string(filename) ) ;
	else
	{
		int index = 0 ;
		float* pbuf = buf ;

		while (index < cube_size) {
			
			if (!fread (pbuf++, sizeof (float), 1, fp))
				throw ReadDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	}
}

void read_my_double_cube_data( const char* filename, int cube_size, double * buf )
{
	FILE * fp = fopen( filename, "r" ) ;
	
	if( !fp ) 
		throw ErrnoError( std::string(filename) ) ;
	else
	{
		int index = 0 ;
		double* pbuf = buf ;

		while (index < cube_size) {
			if (!fread (pbuf++, sizeof (double), 1, fp))
				throw ReadDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	}
}

void write_my_cube_hdr( const char* filehead, int length, int width, int height )
{
	std::string filename = std::string( filehead ) + ".hdr" ;
	FILE * fp = fopen( filename.c_str(), "wt" ) ;
  	
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else 
	{
		fprintf( fp, "%d\n%d\n%d\n", height, width, length ) ;
		fclose( fp ) ;
	}
}

void write_my_byte_cube( const char* filehead, int length, int width, int height, unsigned char * buf )
{
	write_my_cube_hdr( filehead, length, width, height ) ;
  	
	std::string filename = std::string( filehead ) + ".u8" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
  	
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else
	{
		int index = 0 ;
		unsigned char* pbuf = buf ;

		while (index < length * width * height) {

			if (!fwrite (pbuf++, sizeof (unsigned char), 1, fp))
				throw WriteDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	} 
}

void write_my_short_cube( const char* filehead, int length, int width, int height, unsigned short * buf )
{
	write_my_cube_hdr( filehead, length, width, height ) ;
	
	std::string filename = std::string( filehead ) + ".i16" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
  	
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else
	{
		int index = 0 ;
		unsigned short* pbuf = buf ;

		while (index < length * height * width) {
			unsigned short val ;
         	swab ((char*)pbuf++, (char*)&val, sizeof (unsigned short)) ;

			if (!fwrite (&val, sizeof (unsigned short), 1, fp))
				throw WriteDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	}
}

void write_my_single_cube( const char* filehead, int length, int width, int height, float * buf )
{
	write_my_cube_hdr( filehead, length, width, height ) ; 
	
	std::string filename = std::string( filehead ) + ".f32" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
 
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else
	{
		int index = 0 ;
		float* pbuf = buf ;

		while (index < length * width * height) {

			if (!fwrite (pbuf++, sizeof (float), 1, fp))
				throw WriteDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	}
}

void write_my_double_cube( const char* filehead, int length, int width, int height, double * buf )
{
	write_my_cube_hdr( filehead, length, width, height ) ; 
	
	std::string filename = std::string( filehead ) + ".f64" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
  	
	if( !fp ) 
		throw ErrnoError( filename ) ;
	else
	{
		int index = 0 ;
		double* pbuf = buf ;

		while (index < length * width * height) {

			if (!fwrite (pbuf++, sizeof (double), 1, fp))
				throw WriteDataError (std::string(filename)) ;
				
			index ++ ;
		}

		fclose( fp ) ;
	}
}
