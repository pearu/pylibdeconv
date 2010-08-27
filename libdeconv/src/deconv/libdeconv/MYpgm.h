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
 * Filename:  MYpgm.h
 */
 

#ifndef MYPGM_H
#define MYPGM_H


#include "MYerror.h"


/***********************************************************************************************
Support Read/Write PGM Image Format:

P5
#comment line
col row
greylevel
binary data array ( unsigned char: 0 < greylevel < 256, short: 255 < greylevel < 65536 ) 

P2 
#comment line
col row
greylevel
ASCII decimal numbers ( unsigned char: 0 < greylevel < 256, short: 255 < greylevel < 65536 )
************************************************************************************************/


class MagicError : public Error
{
	public:
	MagicError( const char* s, const char* magic )
	{
		_error << "Invalid PGM magic in " << s << " : " << magic << "\n" ;
	}
} ;


class HeaderError : public Error
{
	public:
	HeaderError( const char* s, int w, int h, int g )
	{
		_error << "Invalid PGM parameters in " << s << " : width = " << w 
		       << " , height = " << h << " , greylevel = " << g << "\n" ;
	}
} ;


void read_my_pgm_header     ( const char* filename, int& width, int& height, int& greylevel ) ;
void read_my_byte_pgm_data  ( const char* filename, int image_size, unsigned char  * buf ) ;
void read_my_short_pgm_data ( const char* filename, int image_size, unsigned short * buf ) ;


void write_my_byte_pgm   ( const char* filehead, int width, int height, int greylevel, unsigned char  * buf ) ;
void write_my_short_pgm  ( const char* filehead, int width, int height, int greylevel, unsigned short * buf ) ;


#endif   /* include MYpgm.h */
