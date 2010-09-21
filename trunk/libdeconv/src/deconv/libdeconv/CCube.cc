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
 * Filename:  CCube.cc
 */
 

#include "MYcube.h"
#include "CCube.h"

template<> int CCube< unsigned char >::read( const char * filename )
{
	int length, width, height, cube_size ;
	std::string s = std::string( filename ) ;
	std::string::size_type dot_pos = s.find_last_of( "." ) ;
	std::string filehead = s.substr( 0, dot_pos ) ;
	std::string suffix   = s.substr( dot_pos+1 ) ;
	
	if( suffix == std::string( "u8" ) )
	{
		read_my_cube_hdr( filehead, length, width, height ) ;
		cube_size = length * width * height ;	
		init( length, width, height ) ;
		read_my_byte_cube_data( filename, cube_size, _data ) ;
		return 1 ;
	}
	else
	{
		std::cout << "CCube<unsigned char>::read() failed to recognize file suffix : must be u8.\n " ;
		return 0 ;
	}
}



template<> int CCube< unsigned short >::read( const char * filename )
{
	int length, width, height, cube_size ;
	std::string s = std::string( filename ) ;
	std::string::size_type dot_pos = s.find_last_of( "." ) ;
	std::string filehead = s.substr( 0, dot_pos ) ;
	std::string suffix   = s.substr( dot_pos+1 ) ;

	if( suffix == std::string( "u8" ) || suffix == std::string( "i16" ) )
	{
		read_my_cube_hdr( filehead, length, width, height ) ;
		cube_size = length * width * height ;
		init( length, width, height ) ;	
		if( suffix == std::string( "u8" ) )
		{
			unsigned char* buf = (unsigned char*)fftw_malloc (cube_size * sizeof(unsigned char)) ;
//			ByteArray buf( new unsigned char[ cube_size ] ) ;
			read_my_byte_cube_data( filename, cube_size, buf) ;
			for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (unsigned short) buf[i] ;
			fftw_free (buf) ;
			return 1 ;
		}
		else
		{
			read_my_short_cube_data( filename, cube_size, _data ) ;
			return 2 ;
		}
	}
	else
	{
		std::cout << "CCube<unsigned short>::read() failed to recognize file suffix : must be u8 or i16.\n" ;
		return 0 ;
	}
}



template<> int CCube< int >::read( const char * filename )
{
	int length, width, height, cube_size ;
	std::string s = std::string( filename ) ;
	std::string::size_type dot_pos = s.find_last_of( "." ) ;
	std::string filehead = s.substr( 0, dot_pos ) ;
	std::string suffix   = s.substr( dot_pos+1 ) ;
	
	if( suffix == std::string( "u8" ) || suffix == std::string( "i16" ) )
	{
		read_my_cube_hdr( filehead, length, width, height ) ;
		cube_size = length * width * height ;
		init( length, width, height ) ;		
		if( suffix == std::string( "u8" ) )
		{
			unsigned char* buf = (unsigned char*)fftw_malloc (cube_size * sizeof (unsigned char)) ;
//			ByteArray buf( new unsigned char[ cube_size ] ) ;
			read_my_byte_cube_data( filename, cube_size, buf) ;
			for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (int) buf[i] ;
			fftw_free (buf) ;
			return 1 ;
		}
		else
		{
			unsigned short* buf = (unsigned short*)fftw_malloc (cube_size * sizeof (unsigned short)) ;
//			ShortArray buf( new unsigned short[ cube_size ] ) ;
			read_my_short_cube_data( filename, cube_size, buf ) ;
			for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (int) buf[i] ;
			fftw_free (buf) ;
			return 2 ;
		}
	}
	else
	{
		std::cout << "CCube<int>::read() failed to recognize file suffix : must u8 or i16.\n" ;
		return 0 ;
	}
}


	
template<> int CCube< float >::read( const char * filename )
{
	int length, width, height, cube_size ;
	std::string s = std::string( filename ) ;
	std::string::size_type dot_pos = s.find_last_of( "." ) ;
	std::string filehead = s.substr( 0, dot_pos ) ;
	std::string suffix   = s.substr( dot_pos+1 ) ;

	if( suffix == std::string( "u8" ) || suffix == std::string( "i16" ) || suffix == std::string( "f32" ) )
	{
		read_my_cube_hdr( filehead, length, width, height ) ;
		cube_size = length * width * height ;	
		init( length, width, height ) ;	
		if( suffix == std::string( "u8" ) || suffix == std::string( "i16" ) )
		{
			if( suffix == std::string( "u8" ) )
			{
				unsigned char* buf = (unsigned char*)fftw_malloc (cube_size * sizeof (unsigned short)) ;
//				ByteArray buf( new unsigned char[ cube_size ] ) ;
				read_my_byte_cube_data( filename, cube_size, buf ) ;
				for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (float) buf[i] ;
				fftw_free (buf) ;
				return 1 ;
			}
			else
			{
				unsigned short* buf = (unsigned short*)fftw_malloc (cube_size * sizeof (unsigned short)) ;				
//				ShortArray buf( new unsigned short[ cube_size ] ) ;
				read_my_short_cube_data( filename, cube_size, buf) ;
				for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (float) buf[i] ;
				fftw_free (buf) ;
				return 2 ;
			}
		}
		else
		{
			read_my_single_cube_data( filename, cube_size, _data ) ;
			return 3 ;
		}
	}
	else
	{
		std::cout << "CCube<float>::read() failed to recognize file suffix : must be u8 or i16 or f32.\n" ;
		return 0 ;
	}
}



template<> int CCube< double >::read( const char * filename )
{
	int length, width, height, cube_size ;
	std::string s = std::string( filename ) ;
	std::string::size_type dot_pos = s.find_last_of( "." ) ;
	std::string filehead = s.substr( 0, dot_pos ) ;
	std::string suffix   = s.substr( dot_pos+1 ) ;
	
	if( suffix == std::string( "u8" )  || suffix == std::string( "i16" ) || 
	    suffix == std::string( "f32" ) || suffix == std::string( "f64" ) )
	{
		read_my_cube_hdr( filehead, length, width, height ) ;
		cube_size = length * width * height ;
		init( length, width, height ) ;
		if( suffix == std::string( "u8" ) || suffix == std::string( "i16" ) || suffix == std::string( "f32" ) )
		{ 		
			if( suffix == std::string( "u8" ) )
			{
				unsigned char* buf = (unsigned char*)fftw_malloc (cube_size * sizeof (unsigned char)) ;
//				ByteArray buf( new unsigned char[ cube_size ] ) ;
				read_my_byte_cube_data( filename, cube_size, buf ) ;
				for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (double) buf[i] ;
				fftw_free (buf) ;
				return 1 ;
			}
			else if( suffix == std::string( "i16" ) )
			{
				unsigned short* buf = (unsigned short*)fftw_malloc (cube_size * sizeof (unsigned short)) ;
//				ShortArray buf( new unsigned short[ cube_size ] ) ;
				read_my_short_cube_data( filename, cube_size, buf ) ;
				for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (double) buf[i] ;
				fftw_free (buf) ;
				return 2 ;
			}
			else
			{
				float* buf = (float*)fftw_malloc (cube_size * sizeof (float)) ;
//				SingleArray buf( new float[ cube_size ] ) ;
				read_my_single_cube_data( filename, cube_size, buf ) ;
				for( int i = 0 ; i < cube_size ; i++ ) _data[i] = (double) buf[i] ;
				fftw_free (buf) ;
				return 3 ;
			}
		}
		else
		{
			read_my_double_cube_data( filename, cube_size, _data ) ;
			return 4 ;
		}
	}
	else
	{
		std::cout << "CCube<double>::read() failed to recognize file suffix : must be u8 or i16 or f32 or f64.\n" ;
		return 0 ;
	}
}



template<> void CCube< unsigned char >::write( const char * filehead )
{
	if( Valid( true ) ) 
	{
		write_my_byte_cube( filehead, _length, _width, _height, _data ) ;
	}
}



template<> void CCube< unsigned short >::write( const char * filehead )
{
	if( Valid( true ) ) 
	{
		write_my_short_cube( filehead, _length, _width, _height, _data ) ;
	}
}



template<> void CCube< float >::write( const char * filehead )
{
	if( Valid( true ) ) 
	{
		write_my_single_cube( filehead, _length, _width, _height, _data ) ;
	}
}



template<> void CCube< double >::write( const char * filehead )
{
	if( Valid( true ) ) 
	{
		write_my_double_cube( filehead, _length, _width, _height, _data ) ;
	}
}
