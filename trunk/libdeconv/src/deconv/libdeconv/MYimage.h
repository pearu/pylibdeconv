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
 * Filename:  MYimage.h
 */


#ifndef MYIMAGE_H
#define MYIMAGE_H


#include "MYpgm.h"


#define read_my_image_header      read_my_pgm_header
#define read_my_byte_image_data   read_my_byte_pgm_data
#define read_my_short_image_data  read_my_short_pgm_data
#define write_my_byte_image       write_my_byte_pgm
#define write_my_short_image      write_my_short_pgm


#endif   /* include MYimage.h */
