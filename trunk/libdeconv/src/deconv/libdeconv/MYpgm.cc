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
 * Filename:  MYpgm.cc
 */


#include <stdio.h>
#include "MYpgm.h"


void read_my_pgm_header( const char * filename, int & width, int & height, int & greylevel )
{
	char   header[1024], magic[3] ;
	FILE * fp = fopen( filename, "rb" ) ;
	
	if( !fp )
	{
		throw ErrnoError( std::string(filename) ) ;
	}
	else
	{
		fgets( header, 1024, fp ) ;
		magic[0] = header[0] ;
		magic[1] = header[1] ;
		magic[2] = '\0' ;
		if( magic[0] != 'P' || (magic[1] != '5' && magic[1] != '2') )
		{
			throw MagicError( filename, magic ) ;
		}
		else
		{
			fgets( header, 1024, fp ) ;
			while( header[0] == '#' || header[0] == '\n' || header[0] == '\r' || header[0] == '\0' ) 
			{
				fgets( header, 1024, fp ) ;
			}						
			if( sscanf( header, "%d %d", &width, &height ) != 2 )
			{
				width     = 0 ;
				height	  = 0 ;
				greylevel = 0 ;
			}
			else
			{
				fgets( header, 1024, fp ) ;		
				if( sscanf( header, "%d", &greylevel ) != 1 )
				{
					greylevel = 0 ;
				}
			}
		}
		fclose( fp ) ;
		if( width <= 0 || height <= 0 || greylevel < 1 || greylevel > 65535 )
		{
			throw HeaderError( filename, width, height, greylevel ) ;
		}
	}
}



void read_my_byte_pgm_data( const char* filename, int image_size, unsigned char * buf )
{
	char   header[1024], magic ;
	FILE * fp = fopen( filename, "rb" ) ;
	 
	fgets( header, 1024, fp ) ;
	magic = header[1];
	fgets( header, 1024, fp ) ;        
	while( header[0] == '#' || header[0] == '\n' || header[0] == '\r' || header[0] == '\0' )
	{
		fgets( header, 1024, fp ) ;
	}
	fgets( header, 1024, fp ) ;
	if( magic == '5' ) 
	{
		if( ((int) fread( buf, sizeof(unsigned char), image_size, fp )) != image_size )
		{
			throw ReadDataError( std::string(filename) ) ;
		}
	}
	if( magic == '2' )
	{
		int count = 0 ;
		while( !feof( fp ) ) 
		{
			do 
			{
				magic = getc( fp ) ;
			}	while( (magic == ' ' || magic == '\t' || magic == '\n' || magic == '\r') 
				       && !feof( fp ) ) ;
			if( !feof( fp ) )
			{
				unsigned int intensity = 0 ;
				do 
				{
        				unsigned int const digitVal = magic - '0' ;
        				intensity = intensity * 10 + digitVal ;
        				magic = getc( fp ) ;
				}	while( magic >= '0' && magic <= '9' ) ; 
				if( intensity > 255 )
				{
					throw ReadDataError( std::string(filename) ) ;
				}
				else
				{
					buf[count] = (unsigned char) intensity ;
					count++ ;
				}
			}
		}
		if( count != image_size )
		{
			throw ReadDataError( std::string(filename) ) ;
		}
	} 
	fclose( fp ) ;
}



void read_my_short_pgm_data( const char* filename, int image_size, unsigned short * buf )
{
	char   header[1024], magic ;
	FILE * fp = fopen( filename, "rb" ) ;
	 
	fgets( header, 1024, fp ) ;
	magic = header[1];
	fgets( header, 1024, fp ) ;        
	while( header[0] == '#' || header[0] == '\n' || header[0] == '\r' || header[0] == '\0' )
	{
		fgets( header, 1024, fp ) ;
	}
	fgets( header, 1024, fp ) ;
	if( magic == '5' )
	{
		if( ((int) fread( buf, sizeof(unsigned short), image_size, fp )) != image_size )
		{
			throw ReadDataError( std::string(filename) ) ;
		}
	}
	if( magic == '2' )
	{
		int count = 0 ;
		while( !feof( fp ) ) 
		{
			do 
			{
				magic = getc( fp ) ;
			}	while( (magic == ' ' || magic == '\t' || magic == '\n' || magic == '\r') 
				       && !feof( fp ) ) ;
			if( !feof( fp ) )
			{
				unsigned int intensity = 0 ;
				do 
				{
        				unsigned int const digitVal = magic - '0' ;
        				intensity = intensity * 10 + digitVal ;
        				magic = getc( fp ) ;
				}	while( magic >= '0' && magic <= '9' ) ; 
				if( intensity > 65535 )
				{
					throw ReadDataError( std::string(filename) ) ;
				}
				else
				{
					buf[count] = (unsigned char) intensity ;
					count++ ;
				}
			}
		}
		if( count != image_size )
		{
			throw ReadDataError( std::string(filename) ) ;
		}
	}	
	fclose( fp ) ;
}



void write_my_byte_pgm( const char* filehead, int width, int height, int greylevel, unsigned char * buf )
{
	std::string filename = std::string( filehead ) + ".pgm" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
    	
	if( !fp ) 
	{
		throw ErrnoError( filename ) ;
	}
	else 
	{
		fprintf( fp, "P5\n%d %d\n%d\n", width, height, greylevel ) ;
		int image_size = width * height ;
		if( ((int) fwrite( buf, sizeof(unsigned char), image_size, fp )) != image_size )
		{
			throw WriteDataError( filename ) ;
		}
		fclose( fp ) ;
	}
}
		

	
void write_my_short_pgm( const char* filehead, int width, int height, int greylevel, unsigned short * buf )
{
	std::string filename = std::string( filehead ) + ".pgm" ;
	FILE * fp = fopen( filename.c_str(), "wb" ) ;
    	
	if( !fp ) 
	{
		throw ErrnoError( filename ) ;
	}
	else 
	{	
		fprintf( fp, "P5\n%d %d\n%d\n", width, height, greylevel ) ;
		int image_size = width * height ;		
		if( ((int) fwrite( buf, sizeof(unsigned short), image_size, fp )) != image_size )
		{
			throw WriteDataError( filename ) ;
		}
		fclose( fp ) ;
	}
}
