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
 * Filename:  SHIFTfft.h
 */
 

#ifndef SHIFTFFT_H
#define SHIFTFFT_H


void shift1d( int dimX, float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;
void shift1d( int dimX, double * in_re, double * in_im, double * out_re, double * out_im ) ;


void shift2d( int dimX, int dimY, float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;
void shift2d( int dimX, int dimY, double * in_re, double * in_im, double * out_re, double * out_im ) ;


void shift3d( int dimX, int dimY, int dimZ, float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;
void shift3d( int dimX, int dimY, int dimZ, double * in_re, double * in_im, double * out_re, double * out_im ) ;


void shift( int rank, int dimX, int dimY, int dimZ, float  * in_re, float  * in_im, float  * out_re, float  * out_im ) ;
void shift( int rank, int dimX, int dimY, int dimZ, double * in_re, double * in_im, double * out_re, double * out_im ) ;


#endif /*   include  "SHIFTfft.h"  */
