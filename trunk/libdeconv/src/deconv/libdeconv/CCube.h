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
 * Filename:  CCube.h
 */


#ifndef CCUBE_H_DEFINED
#define CCUBE_H_DEFINED


#include "CSlice.h"


/* 
 *	The CCube template class is a cubic image container designed for handling 3-D greylevel images.
 *
 *      Length (X) is the fastest varying dimension of a cubic image and it indexes from x = 0 to length()-1.
 *	Width  (Y) is the middle          dimension of a cubic image and it indexes from y = 0 to width()-1.
 *	Height (Z) is the slowest varying dimension of a cubic image and it indexes from z = 0 to height()-1.
 *
 *	The cubic image data is stored in the CCube one-dimensional <DataArray> as "x + y*Length()+ z*Width()*Length()".
 *
 *
 *	The CCube can be saved in two files: 	
 *						a header file whose suffix is "hdr" and 
 *						a data file whose suffix depends on the storing data.
 *	The header file is to hold the CCube dimensions and written in the following format: 
 *						<The height of the CCube> //the slowest varying dimension
 *						<the width  of the CCube>
 *						<the length of the CCube> //the fast varying dimension 
 *	The data file is to store the CCube one-dimensional <DataArray> in a binary stream:
 *						Its suffix is "u8"  if data is 8-bit  unsigned char.
 *						Its suffix of "i16" if data is 16-bit unsigned int. 
 *						Its suffix of "f32" if data is 32-bit floating.	
 *						Its suffix of "f64" if data is 64-bit double.
 */


#define XY_SLICE 1	// Indicator of fetching 2-D slices along the Z-dimension or height from the CCube. 
#define XZ_SLICE 2	// Indicator of fetching 2-D slices along the Y-dimension or width  from the CCube.
#define YZ_SLICE 3	// Indicator of fetching 2-D slices along the X-dimension or length from the CCube.



class CubeError : public Error
{
	public:
	CubeError( int length, int width, int height, void* ptr )
	{
		_error << "Invalid cube : ("
		       << length
		       << "x"
		       << width
		       << "x"
		       << height
		       << ") Data : "
		       << ptr
		       << "\n" ;
	}
} ;


class CubeSizeError : public Error
{
	public:
	CubeSizeError( int length, int width, int height )
	{
		_error << "Invalid cube dimensions : "
		       << length
		       << "x"
		       << width
		       << "x"
		       << height
		       << "\n" ;
        }
} ;

	
class CubeAccessError : public Error
{
	public:
	CubeAccessError( int x, int y, int z, int length, int width, int height )
	{
		_error << "Invalid accessing voxel at location : ("
		       << x
		       << ","
		       << y
		       << ","
		       << z
		       << ") in cube size : "
		       << length
		       << "x"
		       << width
		       << "x"
		       << height
		       << "\n" ;
	}
} ;


template < typename T >
class CCube {
	
	
	public:
//	typedef boost::shared_array< T >  DataArray ;
	
	
	/*
	 *	Copy constructor : copy cube to the CCube
	 */
	CCube( const CCube& cube )
	{
		*this = cube ;
	}


	/*
	 *	Constructor
	 *	Input:
	 *		length, it is the length of the CCube; its default value is 0.  
	 *		width,  it is the width  of the CCube; its default value is 0.
	 *		height, it is the height of the CCube; its default value is 0.
	 *		data,   it is an one-dimensional array which stores the CCube data; its default value is NULL.
	 *	Throw:
	 *		throw an error if length < 0 or width < 0 or height < 0.
	 */
	CCube( int length = 0, int width = 0, int height = 0, T* data = NULL )
	{ 
		init( length, width, height, data ) ;
	}
  	

	/*
	 *	Operator "=" constructor : assign cube to the CCube
	 */
	CCube&	operator=( const CCube& cube )
	{
		_length = cube.length() ;
		_width  = cube.width()  ;
		_height = cube.height() ;
//		if( _data ) _data.reset() ;

		if (_data)
			fftw_free(_data) ;
			
//		DataArray temp( new T[ _length * _width * _height ] ) ;
		_data = (T*)fftw_malloc (sizeof(T) *_length * _width * _height);

		for( int i = 0 ; i < _length * _width * _height ; i++ )
		{
			_data[i] = (cube.data())[i] ;
		}
		return *this ;
	}


	/* 
	 *	Operator "==" : compare cube with the CCube
	 *	Output:
	 *		return true  if they are identical.
	 *		return false if they are not identical.
	 */
	bool	operator==( const CCube& cube ) const
	{
		if( _length != cube.length() || _width != cube.width() || _height != cube.height() )       
		{
			return false ;
		}
		else
		{
			for( int i = 0 ; i < _length * _width * _height ; i++ )
			{ 
				if( _data[i] != (cube.data())[i] ) return false ;
			}
			return true ; 
		}
	}  	

	
	/*
	 *	Initialize the CCube 
	 *	Input:
	 *      	length, it is the length of the CCube; its default value is 0.  
	 *		width,  it is the width  of the CCube; its default value is 0.
	 *		height, it is the height of the CCube; its default value is 0.
	 *		data,   it is an one-dimensional array which stores the CCube data; its default value is NULL.
	 *	Throw:
	 *		throw an error if length < 0 or width < 0 or height < 0.
	 */
	void	init( int length = 0, int width = 0, int height = 0, T* data = NULL )
        {
		if( length > 0 && width > 0 && height > 0 )
		{
			_length = length ;
			_width  = width ;
			_height = height ;
			
			if( _data )
				fftw_free(_data) ;

			_data = (T*)fftw_malloc (sizeof(T) *_length * _width * _height);
							
			if (data)
			{
//				DataArray temp( new T[ length * width * height ] ) ;
//				_data = temp ;
				for( int i = 0 ; i < length * width * height ; i++ )
				{
					_data[i] = data[i] ;
				}
			}
		}
		else if( length == 0 && width == 0 && height == 0 )
		{
			_length = 0 ;
			_width  = 0 ;
			_height = 0 ;
			_data = NULL ;
		}	
		else 
		{
			throw CubeSizeError( length, width, height ) ;
		}
	}


	/*	
	 *	Check whether the CCube has been assigned with dimensions or not
	 *	Input:
	 *		throw_if_not, it is to indicate whether to throw an error; its default value is false.
	 *	Output: 
	 *		return true  if the CCube has been assigned with dimensions.
	 * 		return false if the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the input is true and the CCube has not been assigned with dimensions.
	 */
	bool	Valid( bool throw_if_not = false )
	{
		bool ret = ( _length > 0 && _width > 0 && _height > 0 && _data ) ;
                
		if( !ret && throw_if_not ) 
		{
			throw CubeError( _length, _width, _height, _data ) ;
		}
                
		return ret ;
	}

	
	/*	
	 *	Free the CCube : its dimensions will be reset to 0 and data will be freed.
	 */
	void	free()
	{
		if( Valid( false ) )
		{
			_length = 0 ;
			_width  = 0 ;
			_height = 0 ;
			fftw_free(_data) ;
			_data = NULL ;
//			_data.reset() ;
		}
	}

	
	/*
	 *	Get private members
	 *	length() returns the length of the CCube.
	 *	width()  returns the width  of the CCube.
	 *	height() returns the height of the CCube.
	 *	size()   returns the size of the CCube and it equals to length()*width()*height().
	 *	data()   returns the array storing the CCube data.
	 */
  	int	length  () const                   { return _length ;                                   }
  	int	width   () const                   { return _width  ;                                   }
  	int	height  () const                   { return _height ;                                   }
	int	size    () const                   { return ( _length * _width * _height ) ;            }
  	T*	data    () const                   { return _data ;                               }
  	
  	
  	/*
  	 *	Operator "()" : get the value of an voxel in the CCube
  	 *	Input:
  	 *		x, it is the voxel index along the length.
  	 *		y, it is the voxel index along the width.
  	 *		z, it is the voxel index along the height.
  	 *	Output:
  	 *		return the voxel valuce.  
  	 *	Warning:
  	 *		this routine doesnot perform boundary check for each dimension.
  	 */
	T&	operator() ( int x, int y, int z ) 
	{
		return _data[ z*_length*_width + y*_length + x ] ;
	}
	

	/*
  	 *	Set the value of an voxel in the CCube
  	 *	Input:
  	 *		x,     it is the voxel index along the length.
  	 *		y,     it is the voxel index along the width.
  	 *		z,     it is the voxel index along the height.
  	 *		value, it is the voxel value to be set.
  	 *	
  	 *	Throw:
  	 *		throw an error if an index is out of the CCube boundary.
  	 */
	void	setvoxel( int x, int y, int z, T value )
	{
		if( x < 0 || x >= _length || y < 0 || y >= _width || z < 0 || z >= _height )
		{
			throw CubeAccessError( x, y, z, _length, _width, _height ) ;
		}
		else    _data[ z*_length*_width + y*_length + x ] = value ;
	}


	/*
  	 *	Get the value of an voxel in the CCube
  	 *	Input:
  	 *		x,     it is the voxel index along the length.
  	 *		y,     it is the voxel index along the width .
  	 *		z,     it is the voxel index along the height.
  	 *		value, it is to store the voxel value obtained.
  	 *	
  	 *	Throw:
  	 *		throw an error if an index is out of the CCube boundary.
  	 */
	void	getvoxel( int x, int y, int z, T& value )
	{
		if( x < 0 || x >= _length || y < 0 || y >= _width || z < 0 || z >= _height )
		{
			throw CubeAccessError( x, y, z, _length, _width, _height ) ;
		}
		else    value = _data[ z*_length*_width + y*_length + x ] ;
	}


	/*	
	 *	Fill the CCube with a same value
	 *	Input:
	 *		value, it is the value to be filled into the CCube.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	void	fillcube( T value )
	{
		if( Valid( true ) )
		{
			for( int i = 0 ; i < _length * _width * _height ; i++ )
			{
		 		_data[i] = value ;
		 	}
		}
	}


	/*
	 *	Get the data of a 2-D slice from the CCube and put the data into a CSlice
	 *	Input:
	 *		slice_index, it is the index number of the 2-D slice along a CCube dimension.
	 *		slice_plane, it is to indicate along which dimension of the CCube to get the 2-D slice.
	 *		slice,       it is the CSlice to store the 2-D slice which is got from the CCube.
	 *	Output:
	 *		return  0 if success
	 *		return  XY_SLICE(1) if fail for the slice_index is out of the the CCube height.
	 *		return  XZ_SLICE(2) if fail for the slice_index is out of the the CCube width.
	 *		return  YZ_SLICE(3) if fail for the slice_index is out of the the CCube length.
	 *		return  4 if fail for the slice_plane is not one of XY_SLICE(1), XZ_SLICE(2), YZ_SLICE(3).		
	 */
	int	getslice( int slice_index , int slice_plane , CSlice<T>& slice ) ;  
  	
  	
	/*
	 *	Set the data of a slice in the CCube using the data of a 2-D slice stored in a CSlice
	 *	Input:
	 *		slice_index, it is the index number of the 2-D slice along a CCube dimension.
	 *		slice_plane, it is to indicate along which dimension of the CCube to get the 2-D slice.
	 *		slice,       it is the CSlice storing the 2-D slice which will be set into the CCube.
	 *	Output:
	 *		return  0 if success
	 *		return  XY_SLICE (1)  if fail for the slice_index is out of the the CCube height.
	 *		return  XZ_SLICE (2)  if fail for the slice_index is out of the the CCube width.
	 *		return  YZ_SLICE (3)  if fail for the slice_index is out of the the CCube length.
	 *		return  -XY_SLICE(-1) if fail for the CSlice XY-dimensions donot match with CCube XY-dimensions.
	 *		return  -XZ_SLICE(-2) if fail for the CSlice XZ-dimensions donot match with CCube XZ-dimensions.
	 *		return  -YZ_SLICE(-3) if fail for the CSlice YZ-dimensions donot match with CCube YZ-dimensions.
	 *		return  4 if fail for the slice_plane is not one of XY_SLICE(1), XZ_SLICE(2), YZ_SLICE(3).
	 */
	int	setslice( int slice_index , int slice_plane , CSlice<T>  slice ) ;


	/*	
	 *	Draw a homogeneous cyliner in the CCube
	 *	Input:
	 *		cx,    it is the cylinder center along the length of the CCube.
	 *		cy,    it is the cylinder center along the width  of the CCube.
	 *		cz,    it is the cylinder center along the height of the CCube.
	 *		rx,    it is the cylinder radius along the length of the CCube.
	 *		ry,    it is the cylinder radius along the width  of the CCube.
	 *		rz,    it is the cylinder height along the height of the CCube.
	 *		value, it is the cylinder intensity value.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the cylinder is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	draw_cylinder( int cx, int cy, int cz, int rx, int ry, int hz, T value ) ;
	
	
	/*	
	 *	Draw a homogeneous ellipse in the CCube
	 *	Input:
	 *		int cx,  it is the ellipse center along the length of the CCube.
	 *		cy,  it is the ellipse center along the width  of the CCube.
	 *		cz,  it is the ellipse center along the height of the CCube.
	 *		rx,  it is the ellipse radius along the length of the CCube.
	 *		ry,  it is the ellipse radius along the width  of the CCube.
	 *		rz,  it is the ellipse height along the height of the CCube.
	 *		value, it is the ellipse intensity value.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the ellipse is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	draw_ellipse( int cx, int cy, int cz, int rx, int ry, int rz, T value ) ;


	/*
	 *	Read the CCube data from a pair of header (.hdr) and data (.u8/.i16/.f32/.f64) files (described above in this file)
	 *	Input: 
	 *		filename, it is the data file name including the suffix.
	 *	Output:
	 *		return 1 if success in reading a CCube with 8-bit  unsigned   char data.
	 *		return 2 if success in reading a CCube with 16-bit unsigned    int data. 
	 *		return 3 if success in reading a CCube with 32-bit single floating data.
	 *		return 4 if success in reading a CCube with 64-bit double floating data.
	 *		return 0 if fail for the mismatch of the CCube with the CCube data file. 
	 */
	int	read( const char * filename ) ;
  	
  	
	/*
	 *	Write the CCube data to a pair of header(.hdr) and data (.u8/.i16/.f32/.f64) files (described above in this file)
	 *	Input: 
	 *		filehead, it is the file name without the suffix.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	void	write( const char * filehead ) ;


	/*
	 *	Crop the CCube length
	 *	Input:
	 *		index0, it is the lower  length index below which the slices will be cropped.
	 *		index1, it is the higher length index above which the slices will be cropped.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the cropping is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	length_crop( int index0, int index1 ) ;
  	
  	
	/*
	 *	Crop the CCube width
	 *	Input:
	 *		index0, it is the lower  width index below which the slices will be cropped.
	 *		index1, it is the higher width index above which the slices will be cropped.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the cropping is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	width_crop( int index0, int index1 ) ;
  	
  	
	/*
	 *	Crop the CCube height
	 *	Input:
	 *		index0, it is the lower  height index below which the slices will be cropped.
	 *		index1, it is the higher height index above which the slices will be cropped.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the cropping is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	height_crop( int index0, int index1 ) ;
  	
  	
	/*
	 *	Pad the CCube length
	 *	Input:
	 *		num0,  it is the number of slices to be padded before the first slice along the length.
	 *		num1,  it is the number of slices to be padded after  the last  slice along the length.
	 *		value, it is the value assigned to the padded slices.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the padding is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	length_pad( int num0, int num1, T value ) ;
  	
  	
	/*
	 *	Pad the CCube width
	 *	Input:
	 *		num0, it is the number of slices to be padded before the first slice along the width.
	 *		num1, it is the number of slices to be padded after  the last  slice along the width.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the padding is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	width_pad( int num0, int num1, T value ) ;
  	
  	
	/*
	 *	Pad the CCube height
	 *	Input:
	 *		num0,  it is the number of slices to be padded before the first slice along the height.
	 *		num1,  it is the number of slices to be padded after  the last  slice along the height.
	 *		value, it is the value assigned to the padded slices.
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for the padding is out of the CCube boundary.
	 *		return -1 if fail for the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	height_pad( int num0, int num1, T value ) ;


	/*
	 *	Statistical routines
	 *	cubemaxval()  returns the maximum  value of the CCube.
	 *	cubeminval()  returns the minimum  value of the CCube.
	 *	cubemeanval() returns the mean     value of the CCube.
	 *	cubevarval()  returns the variance value of the CCube.
	 *	cubesumval()  returns the summed   value of the CCube.
	 *	throw an error if the CCube has not been assigned with dimensions.
	 */
	double	cubemaxval   () ;
	double	cubeminval   () ;
	double	cubemeanval  () ;
	double	cubevarval   () ;
	double	cubesumval   () ;


	/*
	 *	Shift routine used for shifting the PSF data
	 *	Output:
	 *		return  0 if success.
	 *		return  1 if fail for not all of the CCube dimensions are even.
	 *		return -1 if the CCube has not been assigned with dimensions.
	 *	Throw:
	 *		throw an error if the CCube has not been assigned with dimensions.
	 */
	int	shift() ;   	
  	

	private:
	int        _length ;
	int        _width  ;
	int        _height ;	
	T*  	   _data ;
  	
//	void _init( int length, int width, int height, DataArray data )
	void _init( int length, int width, int height, T* data )
	{
		_length = length ;
		_width  = width ;
		_height = height ;
		_data   = data ;
	}
} ;



template < typename T >
int CCube<T>::getslice( int slice_index , int slice_plane , CSlice<T>& slice )
{
	switch( slice_plane )
	{	
		case XY_SLICE :
			if( slice_index < 0 || slice_index >= _height )
			{ 
				std::cout << "CCbue::getslice() failed for an invalid slice_index : "
				          << slice_index << " to cube-slice indexed from 0 to " << _height << ".\n" ;
				return XY_SLICE ;
			} 
			else 
			{
				slice.init( _length, _width ) ; 
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						slice( i, j ) = _data[ slice_index*_width*_length + j*_length + i ] ;
				return 0 ;
			}
			break ;

		case XZ_SLICE :
			if( slice_index < 0 || slice_index >= _width )
			{
				std::cout << "CCbue::getslice() failed for an invalid slice_index : "
				          << slice_index << " to cube-slice indexed from 0 to " << _width << ".\n" ;
				return XZ_SLICE ;
			}
			else
			{
				slice.init( _length, _height ) ;
				for( int j = 0 ; j < _height ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						slice( i, j ) = _data[ j*_width*_length + slice_index*_length + i ] ;
				return 0 ;
			}
			break ;

		case YZ_SLICE :
			if( slice_index < 0 || slice_index >= _length )
			{
				std::cout << "CCbue::getslice() failed for an invalid slice_index : " 
				          << slice_index << " to cube-slice indexed from 0 to " << _length << ".\n" ;
				return YZ_SLICE ;
			}
			else
			{
				slice.init( _width, _height ) ;
				for( int j = 0 ; j < _height ; ++j )
					for( int i = 0 ; i < _width ; ++i )
						slice( i, j ) = _data[ j*_width*_length + i*_length + slice_index ] ;
				return 0 ;
			}
			break ;
		
		default :
			std::cout << "CCbue::getslice() failed for an invalid slice_plane : "
			          << slice_plane << " (slice_plane must be 1-XY or 2-XZ or 3-YZ).\n" ;
			return 4 ;
			break ;
	}
}



template < typename T >
int CCube<T>::setslice( int slice_index, int slice_plane, CSlice<T> slice )
{
	switch( slice_plane )
	{	
		case XY_SLICE :
			if( slice_index < 0 || slice_index >= _height )
			{ 
				std::cout << "CCbue::setslice() failed for an invalid slice_index : "
				          << slice_index << " to cube-slice indexed from 0 ~ " << _height << ".\n" ;
				return XY_SLICE ;
			} 
			else if( slice.width() != _length || slice.height() != _width )
			{
				std::cout << "CCbue::setslice() failed for an invalid slice : "
				          << slice.width() << "x" << slice.height() 
				          << " to cube-slice( " << _length << "x" << _width << ").\n" ;
				return (-XY_SLICE) ;
			}
			else
			{
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						_data[ slice_index*_width*_length + j*_length + i ] = (slice.data())[ i + j*_length ] ;
				return 0 ;
			}
			break ;

		case XZ_SLICE :
			if( slice_index < 0 || slice_index >= _width )
			{		         
				std::cout << "CCbue::setslice() failed for an invalid slice_index : "
				          << slice_index << " to cube-slice indexed from 0 ~ " << _width << ".\n" ;
				return XZ_SLICE ;
			}
			else if( slice.width() != _length || slice.height() != _height )
			{
				std::cout << "CCbue::setslice() failed for an invalid slice : "
				          << slice.width() << "x" << slice.height() 
				          << " to cube-slice( " << _length << "x" << _height << " ).\n" ;
				return (-XZ_SLICE) ;
			}
			else
			{
				for( int j = 0 ; j < _height ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						_data[ j*_width*_length + slice_index*_length + i ] = (slice.data())[ i + j*_length ] ;
				return 0 ;
			}
			break ;

		case YZ_SLICE :
			if( slice_index < 0 || slice_index >= _length )
			{
				std::cout << "CCbue::setslice() failed for an invalid slice_index : " 
				          << slice_index << " to cube-slice indexed from 0 to " << _length << "\n" ;
				return YZ_SLICE ;
			}
			else if( slice.width() != _width || slice.height() != _height )
			{			 
				std::cout << "CCbue::setslice() failed for an invalid slice : " 
				          << slice.width() << "x" << slice.height() 
				          << " to cube-slice( " << _width << "x" << _height << " ).\n" ;
				return (-YZ_SLICE) ;
			}
			else 
			{
				for( int j = 0 ; j < _height ; ++j )
					for( int i = 0 ; i < _width ; ++i )
						_data[ j*_width*_length + i*_length + slice_index ] = (slice.data())[ i + j*_width ] ;
				return 0 ;
			}
			break ;
    		
		default :		     
			std::cout << "CCbue::getslice() failed for an invalid slice_plane : " 
			          << slice_plane << " (slice_plane must be 1-XY or 2-XZ or 3-YZ).\n" ;
			return 4 ;
			break ;
	}
}



template <typename T>
int CCube<T>::draw_cylinder(int cx, int cy, int cz, int rx, int ry, int hz, T value )
{
	if( Valid( true ) )
	{		  
		if( cx-rx < 0 || cx+rx >= _length || cy-ry < 0 || cy+ry >= _width || 
		    cz-(int)(floor(hz/2)) < 0 || cz+(int)(ceil(hz/2)) >= _height   )
		{
			std::cout << "CCube::draw_cylinder failed for given invalid parameters :\n"
			          << "(cx,cy,cz) = (" << cx << "," << cy << "," << cz << ") , " 
			          << "(rx,ry,hz) = (" << rx << "," << ry << "," << hz << ") "
			          << "to cube (" << _length << "x" << _width << "x" << _height << ").\n" ;
			return 1 ;
		}
		else
		{
			double rx2 = (double) ( rx * rx ) ;
			double ry2 = (double) ( ry * ry ) ;
			for( int k = cz-((int)(floor(hz/2))) ; k <= cz+((int)(ceil(hz/2))) ; ++k )
				for( int j = cy-ry ; j <= cy+ry ; ++j )
					for( int i = cx-rx ; i <= cx+rx ; ++i )
						if( ( ((double)((j-cy)*(j-cy)))/ry2 + ((double)((i-cx)*(i-cx)))/rx2 ) <= 1.0  ) 
							_data[ k*_width*_length + j*_length + i ] = value ;
		}
		return 0 ;
	}
	return -1 ;
}



template <typename T>
int CCube<T>::draw_ellipse(int cx, int cy, int cz, int rx, int ry, int rz, T value )
{
	if( Valid( true ) ) 
	{
		if( cx-rx < 0 || cx-rx >= _length  || cy-ry < 0 || cy-ry >= _width || cz+rz < 0 || cz+rz >= _height )
		{
			std::cout << "CCube::draw_ellips failed for given invalid parameters :\n"
			          << "(cx,cy,cz) = (" << cx << "," << cy << "," << cz << ") , " 
			          << "(rx,ry,rz) = (" << rx << "," << ry << "," << rz << ") "
			          << "to cube (" << _length << "x" << _width << "x" << _height << ").\n" ;
			return 1 ;
		}
		else
		{
			double rx2 = (double) ( rx * rx ) ;
			double ry2 = (double) ( ry * ry ) ;
			double rz2 = (double) ( rz * rz ) ;
			for( int k = cz-rz ; k <= cz+rz ; ++k )
				for( int j = cy-ry ; j <= cy+ry ; ++j )
					for( int i = cx-rx ; i <= cx+rx ; ++i )
						if( ( ((double)((k-cz)*(k-cz)))/rz2 + ((double)((j-cy)*(j-cy)))/ry2 
						                         + ((double)((i-cx)*(i-cx)))/rx2 ) <= 1.0 ) 
							_data[ k*_width*_length + j*_length + i ] = value ;
		}
		return 0 ;
	}
	return -1 ;
}



template <typename T>
int CCube<T>::read( const char * filename ) 
{
	std::cout << filename << " failed to recognize the class data type : Only support\n"
	          << "CCube<unsigned char> | CCube<unsigned short> | CCube<int> | CCube<float> | CCube<double>.\n" ;
	return 0 ;
}
template<> int CCube< unsigned char  >::read( const char * filename ) ;
template<> int CCube< unsigned short >::read( const char * filename ) ;
template<> int CCube< int            >::read( const char * filename ) ;
template<> int CCube< float          >::read( const char * filename ) ;
template<> int CCube< double         >::read( const char * filename ) ;



template <typename T>
void CCube<T>::write( const char* filehead )
{
	std::cout << "CCube::write() to " << filehead << " failed to recognize the class data type : Only support \n"
	          << "CCube<unsigned char> | CCube<unsigned short> | CCube<float> | CCube<double>.\n" ;
}
template<> void CCube< unsigned char  >::write( const char * filehead ) ;
template<> void CCube< unsigned short >::write( const char * filehead ) ;
template<> void CCube< float          >::write( const char * filehead ) ;
template<> void CCube< double         >::write( const char * filehead ) ;



template <typename T>
int CCube<T>::length_crop( int index0, int index1 )
{
	if( Valid( true ) ) 
	{ 
		if( index0 < 0 || index0 > _length-2 || index1 < 1 || index1 > _length-1 || index0 >= index1  )
		{
			std::cout << "CCube::length_crop() failed for given invalid parameters :\n"
			          << "lower_index = " << index0 << " , higher_index = " << index1 
			          << " to cube's length indexed from 0 to " << _length << ".\n" ;
			return 1 ;
		}
		else
		{
			int length = index1 - index0 + 1 ;
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ length * _width * _height ] ) ;			

			for( int k = 0 ; k < _height ; k++ )
				for( int j = 0 ; j < _width ; j++ )
					for( int i = 0 ; i < length ; i++ )
						buf[ k*_width*length + j*length + i ] 
						= _data[ k*_width*_length + j*_length + i+index0 ] ;  	    	 

			_init( length, _width, _height, buf ) ;
			fftw_free (buf) ;

			return 0 ;
		}
	}
	return -1 ;
}

 

template <typename T>   
int CCube<T>::width_crop( int index0, int index1 )
{
	if( Valid( true ) ) 
	{ 
		if( index0 < 0 || index0 > _width-2 || index1 < 1 || index1 > _width-1 || index0 >= index1 )
		{
			std::cout << "CCube::width_crop() failed for given invalid parameters :\n"
			          << "lower_index = " << index0 << " , higher_index = " << index1 
			          << " to cube's width indexed from 0 to " << _width << ".\n" ;
			return 1 ;
		}
		else
		{
			int width = index1 - index0 + 1 ;
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ _length * width * _height ] ) ;

			for( int k = 0 ; k < _height ; k++ )	  	  
				for( int j = 0 ; j < width ; j++ )
					for( int i = 0 ; i < _length ; i++ )
						buf[ k*width*_length + j*_length + i ] 
						= _data[ k*_width*_length + (j+index0)*_length + i ] ;

			_init( _length, width, _height, buf ) ;
			fftw_free (buf) ;

			return 0 ;
		}
	}
	return -1 ;
}



template <typename T>
int CCube<T>::height_crop( int index0, int index1 )
{
	if( Valid( true ) ) 
	{ 
		if( index0 < 0 || index0 > _height-2 || index1 < 1 || index1 > _height-1 || index0 >= index1 )
		{
			std::cout << "CCube::height_crop() failed given invalid parameters :\n"
			          << "lower_index = " << index0 << " , higher_index = " << index1 
			          << " to cube's height indexed from 0 to " << _height << ".\n" ;
			return 1 ;
		}
		else
		{
			int height = index1 - index0 + 1 ;  
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ _length * _width * height ] ) ;

			for( int k = 0 ; k < height ; ++k )	  	  
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						buf[ k*_width*_length + j*_length + i ] 
						= _data[ (k+index0)*_width*_length + j*_length + i ] ;

			_init( _length, _width, height, buf ) ;
			fftw_free (buf) ;

			return 0 ;
		}
	}
	return -1 ;	
}



template <typename T>
int CCube<T>::length_pad( int num0, int num1, T value )
{
	if( Valid( true ) ) 
	{ 
		if( num0 < 0 || num1 < 0 )
		{
			std::cout << "CCube::length_pad() failed given invalid parameters :\n"
			          << "num_slices_add_before = " << num0 << " , num_slices_add_after = " << num1 << ".\n" ;
			return 1 ;
		}
		else
		{
			int length = num1 + num0 + _length ;
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ length * _width * _height ] ) ;
			for( int i = 0 ; i < length * _width * _height ; i++ ) buf[i] = value ;

			for( int k = 0 ; k < _height ; ++k )	  	  
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						buf[ k*_width*length + j*length + i+num0 ] 
						= _data[ k*_width*_length + j*_length + i ] ;

			_init( length, _width, _height, buf ) ;
			fftw_free (buf) ;
			return 0 ;
		}
	}
	return -1 ;
}



template <typename T>
int CCube<T>::width_pad(int num0, int num1, T value)
{
	if( Valid( true ) ) 
	{ 
		if( num0 < 0 || num1 < 0 )
		{
			std::cout << "CCube::width_pad() failed for given invalid parameters :\n"
			          << "num_slices_add_before = " << num0 << " , num_slices_add_after = " << num1 << ".\n" ;
			return 1 ;
		}
		else
		{
			int width = num1 + num0 + _width ;
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ _length * width * _height ] ) ;
			for( int i = 0 ; i < _length * width * _height ; i++ ) buf[i] = value ;

			for( int k = 0 ; k < _height ; ++k )	  	  
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						buf[ k*width*_length + (j+num0)*_length + i ] 
						= _data[ k*_width*_length + j*_length + i ] ;

			_init( _length, width, _height, buf ) ;

			return 0 ;
		}
	}
	return -1 ;
}



template <typename T>
int CCube<T>::height_pad(int num0, int num1, T value)
{
	if( Valid( true ) ) 
	{ 
		if( num0 < 0 || num1 < 0 )
		{
			std::cout << "CCube::height_pad() failed for given invalid parameters :\n"
			          << "num_slices_add_before = " << num0 << " , num_slices_add_after = " << num1 << ".\n" ;
			return 1 ;
		}
		else
		{
			int height = num1 + num0 + _height ;
			T* buf = (T*)fftw_malloc (sizeof(T) * _length * _width * _height) ;
//			DataArray buf( new T [ _length * _width * height ] ) ;
			for( int i = 0 ; i < _length * _width * height ; i++ ) buf[i] = value ;

			for( int k = 0 ; k < _height ; ++k )	  	  
				for( int j = 0 ; j < _width ; ++j )
					for( int i = 0 ; i < _length ; ++i )
						buf[ (k+num0)*_width*_length + j*_length + i ] 
						= _data[ k*_width*_length + j*_length + i ] ;

			_init( _length, _width, height, buf ) ;

			return 0 ;
		}
	}
	return -1 ;
}



template <typename T>
double CCube<T>::cubemaxval () 
{
	if( Valid( true ) ) 
	{
		double max_value = (double) _data[0] ;		

		for( int i = 0 ; i < _length * _width * _height ; ++i )
		{ 
			if( ((double)_data[i]) > max_value ) max_value = (double) _data[i] ;
		}

		return max_value ;
	}
	return -1 ;
}



template <typename T>
double CCube<T>::cubeminval()
{
	if( Valid( true ) ) 
	{
		double min_value = (double) _data[0] ;  		

		for( int i = 0 ; i < _length * _width * _height ; ++i )
		{ 
			if( ((double)_data[i]) < min_value ) min_value = (double) _data[i] ;
		}

		return min_value ;
	}
	return -1 ;
}



template <typename T>
double CCube<T>::cubemeanval() 
{
	if( Valid( true ) ) 
	{
		double mean_value = 0.0 ;  		

		for( int i = 0 ; i < _length * _width * _height ; ++i ) 
		{
			mean_value += ((double)_data[i]) ;
		}
		mean_value /= ((double)(_length * _width * _height)) ;

		return mean_value ;
	}
	return -1 ;
}



template <typename T>
double CCube<T>::cubevarval() 
{
	if( Valid( true ) ) 
	{
		double var_value  = 0.0 ;
		double mean_value = 0.0 ;  		

		for( int i = 0 ; i < _length * _width * _height ; ++i )
		{ 
			mean_value += ((double)_data[i]) ;
		}
		mean_value /= ((double)(_length * _width * _height)) ;	  	
		for( int i = 0 ; i < _length * _width * _height ; ++i )
		{ 
			var_value += ( ( (double)_data[i] - mean_value ) * ( (double)_data[i] - mean_value ) ) ;
		}
		var_value =  sqrt( var_value ) / ( (double)(_length * _width * _height) - 1.0 ) ;  		
	  	
		return var_value ;
	}
	return -1 ;
}



template <typename T>
double CCube<T>::cubesumval() 
{
	if( Valid( true ) ) 
	{
		double sum_value  = 0.0 ;
	  	
		for( int i = 0 ; i< _length * _width * _height ; ++i )
		{ 
			sum_value += ((double)_data[i]) ; 
		}
  		
		return sum_value ;
	}
	return -1 ;
}



template <typename T>
int CCube<T>::shift()
{ 	
	T temp ;
	if( Valid( true ) )
	{
		if( _length%2 == 0 && _width%2 == 0 && _height%2 == 0 ) 
		{
			for( int k = 0 ; k < _height/2 ; k++ ) 
			{
				for( int j = 0 ; j < _width/2 ; j++ ) 
				{
					for( int i = 0 ; i < _length/2 ; i++ ) 
					{
						temp = _data[ i + j*_length + k*_length*_width ] ;
						_data[ i + j*_length + k*_length*_width ] 
						= _data[ i+_length/2 + (j+_width/2)*_length + (k+_height/2)*_length*_width ] ;
						_data[ i+_length/2 + (j+_width/2)*_length + (k+_height/2)*_length*_width ] = temp ;
					}
					for( int i = _length/2 ; i < _length ; i++ ) 
					{
						temp = _data[ i + j*_length + k*_length*_width ] ;
						_data[ i + j*_length + k*_length*_width ] 
						= _data[ i-_length/2 + (j+_width/2)*_length + (k+_height/2)*_length*_width ] ;
						_data[ i-_length/2 + (j+_width/2)*_length + (k+_height/2)*_length*_width ] = temp ;
					}
				}
				for( int j = _width/2 ; j < _width ; j++ ) 
				{
					for( int i = 0 ; i < _length/2 ; i++ ) 
					{
						temp = _data[ i + j*_length + k*_length*_width ] ;
						_data[ i + j*_length + k*_length*_width ] 
						= _data[ i+_length/2 + (j-_width/2)*_length + (k+_height/2)*_length*_width ] ;
						_data[ i+_length/2 + (j-_width/2)*_length + (k+_height/2)*_length*_width ] = temp ;
					}
					for( int i = _length/2 ; i < _length ; i++ ) 
					{
						temp = _data[ i + j*_length + k*_length*_width ] ;
						_data[ i + j*_length + k*_length*_width ] 
						= _data[ i-_length/2 + (j-_width/2)*_length + (k+_height/2)*_length*_width ] ;
						_data[ i-_length/2 + (j-_width/2)*_length + (k+_height/2)*_length*_width ] = temp ;
					}
				}
			}
		}
		else
		{
			std::cout << "CCube::shift() failed for each cube's dimension must be even : "
			          << _length << "x" << _width << "x" << _height << ".\n" ;
			return 1 ;
		}
		return 0 ;
	}
	return -1 ;
}


typedef CCube<unsigned char>  ByteCube ;
typedef CCube<unsigned short> ShortCube ;
typedef CCube<float>          SingleCube ;
typedef CCube<double>         DoubleCube ;


#endif  /* include "CCube.h" */
