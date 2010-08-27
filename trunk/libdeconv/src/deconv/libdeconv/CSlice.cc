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
 * Filename:  CSlice.cc
 */


#include "CSlice.h"



template<> int CSlice< unsigned char >::read( const char * filename )
{
	int width, height, greylevel ;

	read_my_image_header( filename, width, height, greylevel ) ;
	if( greylevel < 256 )
	{
		init( width, height ) ;
		read_my_byte_image_data( filename, width*height, _data.get() ) ;
		return 255 ;
	}
	else
	{
		std::cout << "CSlice::read() failed to read data from " 
		          << filename << " : CSlice<unsigned char> cannot read unsigned short data.\n" ;
		return 0 ;
	}
}



template<> void CSlice< unsigned char >::write( const char* filehead )
{
	if( Valid( true ) )
	{
		 write_my_byte_image( filehead, _width, _height, 255, _data.get() ) ;  
	}
}



template<> void CSlice< unsigned short >::write( const char* filehead )
{
	if( Valid( true ) )
	{
		 write_my_short_image( filehead, _width, _height, 65535, _data.get() ) ;  
	}
}
