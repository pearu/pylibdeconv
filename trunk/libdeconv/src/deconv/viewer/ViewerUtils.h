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
 * Filename:  ViewerUtils.h
 */


#ifndef VIEWERUTILS_H
#define VIEWERUTILS_H


/* 
 *	Check an integer whether it is power of 2 and return true if it is, false otherwise.
 */
bool IsNotPowerOf2( int num ) ;


/* 
 *	Name a sequence of files without suffix
 *	for example : setname = "myfile", ID = 5
 *	              if      EndID < 10   then filehead = "myfile5"
 *	              if  9 < EndID < 100  then filehead = "myfile05"
 *	              if 99 < EndID < 1000 then filehead = "myfile005"
 *	              ... ... ...
 */ 
void name_seq_files( char* filehead, char* setname, int ID, int EndID = 0 ) ;


/*
 *	Write a byte-image buffer to a PPM file
 *	P6                     // magic
 *	width length           // image size
 *	255                    // greylevel
 *	binary stream holding image data, stream size in unit of unsigned char = image size * 3 
 */
void write_my_byte_ppm( char* filehead, int width, int height, unsigned char * buf ) ;


/*
 *	Contrasting functions using the "Window Width" and "Window Center" method 
 */
unsigned char byte_contrast( double intensity, double IntensityWindow, double IntensityCenter ) ;
unsigned short short_contrast( double intensity, double IntensityWindow, double IntensityCenter ) ;
double contrast( double intensity, double IntensityWindow, double IntensityCenter ) ;


/*
 *	Read parameters to construct the shape of a Fourier filter
 */
void set_cylinder( int & rx, int & ry, int & hz ) ;
void set_ellipse ( int & rx, int & ry, int & rz ) ;


#endif /* include ViewerUtils.h */

