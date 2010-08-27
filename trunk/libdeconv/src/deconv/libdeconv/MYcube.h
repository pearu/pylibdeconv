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
 * Filename:  MYcube.h
 */


#ifndef MYCUBE_H
#define MYCUBE_H

#include "MYerror.h"


void read_my_cube_hdr        ( std::string filehead, int& length, int& width, int& height ) ;
void read_my_byte_cube_data  ( const char* filename, int cube_size, unsigned char  * buf ) ;
void read_my_short_cube_data ( const char* filename, int cube_size, unsigned short * buf ) ;
void read_my_single_cube_data( const char* filename, int cube_size, float  * buf ) ;
void read_my_double_cube_data( const char* filename, int cube_size, double * buf ) ;


void write_my_cube_hdr   ( const char* filehead, int length, int width, int height ) ;
void write_my_byte_cube  ( const char* filehead, int length, int width, int height, unsigned char  * buf ) ;
void write_my_short_cube ( const char* filehead, int length, int width, int height, unsigned short * buf ) ;
void write_my_single_cube( const char* filehead, int length, int width, int height, float  * buf ) ;
void write_my_double_cube( const char* filehead, int length, int width, int height, double * buf ) ;


#endif   /* include MYcube.h */
