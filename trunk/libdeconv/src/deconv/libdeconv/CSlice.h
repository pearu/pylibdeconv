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
 * Filename:  CSlice.h
 */


#ifndef CSLICE_H
#define CSLICE_H


#include <math.h>
#include "MYdef.h"
#include "MYimage.h"


/* 
 *	The CSlice template class is an image container designed for handling 2-D greylevel images.
 *
 *	Width  (X) is the fast varying dimension of an image and it indexes from x = 0 to length()-1.
 *	Height (Y) is the slow varying dimension of an image and it indexes from y = 0 to width()-1.
 *
 *	The image data is stored in the CSlice one-dimensional <DataArray> as "x + y * width()".
 *
 *	The CSlice only supports reading from or writing to 8-bit/16-bit PGM images.
 *
 */


class SliceError : public Error
{
	public:	
	SliceError( int width, int height, void* ptr )
	{
		_error << "Invalid slice : ("
		       << width
		       << "x"
		       << height
		       << ") Data : "
		       << ptr
 		       << "\n" ;
	}
} ;
	

class SliceSizeError : public Error
{
	public:
        SliceSizeError( int width, int height )
	{
		_error << "Invalid slice dimensions : "
		       << width
		       << "x"
		       << height
		       << "\n" ;
	}
} ;


class SliceAccessError : public Error
{
	public:
	SliceAccessError( int x, int y, int width, int height )
	{
		_error << "Invalid accessing pixel at location : ("
		       << x
		       << ","
		       << y
		       << ") in slice size : "
		       << width
		       << "x"
		       << height
		       << "\n" ;
	}
} ;


template< typename T >
class CSlice 
{	
	public:	
	typedef boost::shared_array< T >  DataArray ;
	
	
	/*
	 *	Copy constructor : copy slice to the CSlice
	 */
	CSlice( const CSlice& slice )
	{
		*this = slice ;
	}


	/*
	 *	Constructor
	 *	Input:
	 *		width,  it is the width  of the CSlice; its default value is 0.
	 *		height, it is the height of the CSlice; its default value is 0.
	 *		data,   it is an one-dimensional array which stores the CSlice data; its default value is NULL.
	 *	Throw:
	 *		throw an error if width < 0 or height < 0.
	 */
	CSlice( int width = 0, int height = 0, T* data = NULL )
	{
		init( width, height, data) ;
	}

	
	/*
	 *	Operator "=" constructor : assign slice to the CSlice
	 */
	CSlice& operator=( const CSlice& slice )
	{
		_width  = slice.width()  ;
		_height = slice.height() ;
		if( _data ) _data.reset() ;
		
		DataArray temp( new T[ _width * _height ] ) ;
		_data = temp ;
		
		for( int i = 0 ; i < _width * _height ; i++ )
		{
			_data[i] = (slice.data())[i] ;
		}
		return *this ;
	}


	/* 
	 *	Operator "==" : compare slice with the CSlice
	 *	Output: boolean
	 *		return true  if they are identical.
	 *		return false if they are not identical.
	 */
	bool	operator==( const CSlice& slice ) const
	{
		if( _width != slice.width() || _height != slice.height() ) 
		{
			return false ;
		}
		else
		{
			for( int i = 0 ; i < _width * _height ; i++ )
			{ 
				if( _data[i] != (slice.data())[i] ) return false ;
			}
			return true ;
		}
	}


	/*
	 *	Initialize the CSlice
	 *	Input:
	 *		width,  it is the width  of the CSlice; its default value is 0.
	 *		height, it is the height of the CSlice; its default value is 0.
	 *		data,   it is an one-dimensional array which stores the CSlice data; its default value is NULL.
	 *	Throw:
	 *		throw an error if width < 0 or height < 0.
	 */
	void	init( int width = 0, int height = 0,  T* data = NULL )
	{
		if( width > 0 && height > 0 )
		{
			_width  = width ;
			_height = height ;
			if( _data ) _data.reset() ;
			
			if( data == NULL )
			{
				DataArray temp( new T[ width * height ] ) ;
				_data = temp ;
			}
			else 
			{
				DataArray temp( new T[ width * height ] ) ;
				_data = temp ;					
				for( int i = 0 ; i < width * height ; i++ )
				{
					_data[i] = data[i] ;
				}					
			}
		}
		else if( width == 0 && height == 0 )
		{
			_width  = 0 ;
			_height = 0 ;
		}			
		else 
		{
			throw SliceSizeError( width, height ) ;
		}
	}	


	/*	
	 *	Check whether the CSlice has been assigned with dimensions or not
	 *	Input:
	 *		throw_if_not, it is to indicate whether to throw an error; its default value is false.
	 *	Output:
	 *		return true  if the CSlice has been assigned with dimensions.
	 * 		return false if the CSlice has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the input is true and the CSlice has not been assigned with dimensions.
	 */
	bool	Valid( bool throw_if_not = false )
	{
		bool ret = ( _width > 0 && _height > 0 && _data ) ;
		
		if( !ret && throw_if_not ) 
		{
			throw SliceError( _width, _height, _data.get() ) ;
		}
		
		return ret ;
	}


	/*	
	 *	Free the CSlice : its dimensions will be reset to 0 and data will be freed.
	 */
	void	free()
	{
		if( Valid( false ) )
		{
			_width  = 0 ;
			_height = 0 ;
			_data.reset() ;
		}
	}


	/*
	 *	Get private members
	 *	width()  returns the width  of the CSlice.
	 *	height() returns the height of the CSlice.
	 *	size()   returns the size of the CSlice and it equals to width()*height().
	 *	data()   returns the array storing the CSlice data.
	 */
	int	width   () const            { return _width ;               }
	int	height  () const            { return _height ;              }
	int	size    () const            { return ( _width * _height ) ; }
	T*	data    () const            { return _data.get() ;          }
        
        
	/*
	 *	Operator "()" : get the value of an pixel in the CSlice
	 *	Input:
	 *		x, it is the pixel index along the width.
	 *		y, it is the pixel index along the height.
	 *	Output: T
	 *		return the pixel valuce.  
	 *	Warning:
	 *		this routine doesnot perform boundary check for each dimension.
	 */
        T&	operator() ( int x, int y ) 
	{ 
		return _data[ y*_width + x ] ;    
	}


	/*
	 *	Set the value of an pixel in the CSlice
	 *	Input:
	 *		x,     it is the pixel index along the width.
	 *		y,     it is the pixel index along the height.
	 *		value, it is the pixel value to be set.
	 *	Throw:
	 *		throw an error if an index is out of the CCube boundary.
	 */
	void	setpixel( int x, int y, T value )
	{
		if( x < 0 || x >= _width || y < 0 || y >= _height ) 
		{
			throw SliceAccessError( x, y, _width, _height ) ;
		}
		else    _data[ y*_width + x ] = value ;
	}


	/*
	 *	Get the value of an pixel in the CSlice
	 *	Input:
	 *		x,     it is the pixel index along the width.
	 *		y,     it is the pixel index along the height.
	 *		value, it is to store the pixel value obtained.
	 *	Throw:
	 *		throw an error if an index is out of the CCube boundary.
	 */
	void	getpixel( int x, int y, T& value )
	{
		if( x < 0 || x >= _width || y < 0 || y >= _height ) 
		{
			throw SliceAccessError( x, y, _width, _height ) ;
		}
		else    value = _data[ y*_width + x ] ;
	} 


	/*	
	 *	Fill the CSlice with a same value
	 *	Input:
	 *		value, it is the value to be filled into the CSlice.
	 *	Throw:
	 *		throw an error if the CSlice has not been assigned with dimensions.
	 */
	void	fillslice( T value )
	{
		if( Valid( true ) ) 
		{
			for( int i = 0 ; i < _width * _height ; i++ )
			{
				_data[i] = value ;
			}
		}
	}


	/*
	 *	Read the CSlice from a PGM file
	 *	Input: 
	 *		filename, it is the file name including the suffix.
	 *	Output:
	 *		return 255   if success in reading a 8-bit  PGM image.
	 *		return 65535 if success in reading a 16-bit PGM image.
	 *		return 0     if fail
	 */
	int	read( const char * filename ) ;
	
	
	/*
	 *	Write the CSlice to a PGM file
	 *	Input: 
	 *		filehead, it is the file name without the suffix.
	 *	Throw:
	 *		throw an error if the CSlice has not been assigned with dimensions.
	 */
	void	write( const char * filehead ) ;


	/*
	 *	Crop the CSlice to the size newx*newy with respect to CSlice(cx,cy) 
	 *	Input:
	 *		newx,  it will be the width  of the cropped CSlice.
	 *		newy,  it will be the height of the cropped CSlice.
	 *		cx,cy, it will be the center of the cropped CSlice. 
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the cropping is out of the CSlice boundary.
	 *		return -1 if fail for the CSlice has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	crop( int newx, int newy, int cx, int cy ) ;
	
	
	/*
	 *	Pad the CSlice to the size newx*newy with respect to the CSlice center 
	 *	Input:
	 *		newx,  it will be the width  of the padded CSlice.
	 *		newy,  it will be the height of the padded CSlice.
	 *		value, it will be the value  of the padded pixels.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the padding is out of the CSlice boundary.
	 *		return -1 if fail for the CSlice has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	pad( int newx, int newy, T value ) ;


	/*
	 *	Statistical routines
	 *	slicemaxval()  returns the maximum  value of the CSlice.
	 *	sliceminval()  returns the minimum  value of the CSlice.
	 *	slicemeanval() returns the mean     value of the CSlice.
	 *	slicevarval()  returns the variance value of the CSlice.
	 *	slicesumval()  returns the summed   value of the CSlice.
	 *	throw an error if the CCube has not been assigned with dimensions.
	 */
	double slicemaxval  () ;
	double sliceminval  () ;
	double slicemeanval () ;
	double slicevarval  () ;
	double slicesumval  () ; 


	private :		
	int       _width ;
	int       _height ;
	DataArray _data ;
	
	void _init( int width, int height, DataArray data )
	{
		_width  = width ;
		_height = height ;
		_data   = data ;
	}
} ;



template <typename T>
int CSlice<T>::read( const char * filename )
{
	int width, height, greylevel ;
	
	read_my_image_header( filename, width, height, greylevel ) ;
	init( width, height ) ;
	if( greylevel < 256 ) 
	{
		ByteArray buf( new unsigned char[ width * height ] ) ;
		read_my_byte_image_data( filename, width*height, buf.get() ) ;
		for( int i = 0 ; i < width*height ; i++ ) _data[i] = (T) buf[i] ;
		return 255 ;
	}
	else 
	{
		ShortArray buf( new unsigned short[ width * height ] ) ;
		read_my_short_image_data( filename, width*height, buf.get() ) ;
		for( int i = 0 ; i < width*height ; i++ ) _data[i] = (T) buf[i] ;
		return 65535 ;
	}
}
template<> int CSlice< unsigned char >::read( const char * filename ) ;



template <typename T>
void CSlice<T>::write( const char * filehead )
{
	std::cout << filehead << " failed to recognize the class data type : Only support\n"
	          << "CSlice<unsigned char> | CSlice<unsigned short>.\n" ;
}
template<> void CSlice< unsigned char  >::write( const char * filehead ) ;
template<> void CSlice< unsigned short >::write( const char * filehead ) ;



template <typename T>
int CSlice<T>::crop( int width, int height, int cx, int cy )
{	
	int half_width, odd_width, half_height, odd_height ;
	
	if( Valid( true ) )
	{
		if( width%2 != 0 )
		{ 
			half_width = (width-1) / 2 ;
			odd_width  = 1 ;
		}
		else
		{
			half_width = width / 2 ;
			odd_width  = 0 ;
		}
		if( height%2 != 0 )
		{
			half_height = (height-1) / 2 ;
			odd_height  = 1 ;
		}
		else
		{
			half_height = height / 2 ;
			odd_height  = 0 ;
		}
		
		if( cx-half_width >= 0 && cx+half_width+odd_width <= _width && 
		    cy-half_height >= 0 && cy+half_height+odd_height <= _height )
		{
			int k = 0 ;
			DataArray temp( new T[ width * height ] ) ;
			
			for( int j = cy-half_height ; j < cy+half_height+odd_height ; j++ )
				for( int i = cx-half_width ; i < cx+half_width+odd_width ; i++ )
				{
					temp[k] = _data[ i + j*_width ] ;
					k++ ;
				}

			_init( width, height, temp ) ;
		}
		else
		{
			std::cout << "CSlice::crop failed to crop " << width << "x" << height 
			          << " in slice(" << _width << "x" << _height << ") at the center (" << cx << "," << cy << ").\n" ;
			return 1 ;
		}
		return 0 ;
	}
	return -1 ;
}



template <typename T>
int CSlice<T>::pad( int width, int height, T value )
{
	if( Valid( true ) )
	{
		if( width >= _width || height >= _height )
		{
			DataArray temp( new T[ width * height ] ) ;
			for( int i = 0 ; i < width * height ; i++ ) _data[i] = value ;
			
			for( int j = 0 ; j < _height ; j++ )
				for( int i = 0 ; i < _width ; i++ )
				{
					temp[ ( j + ((int)(ceil(height/2))) - ((int)(ceil(_height/2))) ) * width +
					      i + ((int)(ceil(width/2))) - ((int)(ceil(_width/2))) ] = _data[ j*_width + i ] ;
				}

			_init( width, height, temp ) ;
		}
		else
		{
			std::cout << "CSlice::pad failed to pad slice(" << _width << "x" << _height << ") to "
			          << width << "x" << height << ".\n" ;
			return 1 ;
		}
		return 0 ;
	}
	return -1 ;
}



template <typename T>
double CSlice<T>::slicemaxval()
{
	if( Valid( true ) )
	{
		double max = (double) _data[0] ;
		
		for( int i = 0 ; i < _width * _height ; i++ )
		{ 
			if( ((double)_data[i]) > max ) max = (double) _data[i] ;
		}	
	        
		return max ;
	}
	return -1 ;
}



template <typename T>
double CSlice<T>::sliceminval()
{
	if( Valid( true ) )
	{
		double min = (double) _data[0] ;
	        
		for( int i = 0 ; i < _width * _height ; i++ )
		{ 
			if( ((double)_data[i]) < min ) min = (double) _data[i] ;
		}
		
		return min ;
	}
	return -1 ;
}



template <typename T>
double CSlice<T>::slicemeanval()
{
	if( Valid( true ) )
	{	
		double mean = 0.0 ;
		
		for( int i = 0 ; i < _width * _height ; i++ )
		{
			mean += ( (double) _data[i] ) ;
		}
		mean /= ( (double) (_width * _height) ) ;
	        
		return mean ;
	}
	return -1 ;
}



template <typename T>
double CSlice<T>::slicevarval()
{
	if( Valid( true ) )
	{
		double mean = 0.0 ;
		double var  = 0.0 ;
		
		for( int i = 0 ; i < _width * _height ; i++ )
		{ 
			mean += ( (double) _data[i] ) ;
		}
		mean /= ( (double) (_width * _height) ) ;
		for( int i = 0 ; i < _width * _height ; i++ )
		{ 
			var += ( ( (double) _data[i] - mean ) * ( (double) _data[i] - mean ) ) ;
		}
		var =  sqrt(var) / ( (double) (_width * _height) - 1.0 ) ;
		
		return var ;
	}
	return -1 ;
}



template <typename T>
double CSlice<T>::slicesumval()
{
	if( Valid( true ) )
	{
		double sum = 0.0 ;
		
		for( int i = 0 ; i < _width * _height ; i++ ) 
		{
			sum += ( (double) _data[i] ) ;
		}
	        
		return sum ;
	}
	return -1 ;
}


typedef CSlice<unsigned char>  ByteSlice ;
typedef CSlice<unsigned short> ShortSlice ;
typedef CSlice<float>          SingleSlice ;
typedef CSlice<double>         DoubleSlice ;

	
#endif  /* include CSlice.h */
