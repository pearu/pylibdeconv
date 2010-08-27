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


#include <stdlib.h>
#include <math.h>
#include "SHIFTfft.h"


void shift1d( int dimX, float * in_re, float * in_im, float * out_re, float * out_im )
{
	if( in_im == NULL )
	{
		for( int i = 0 ; i < dimX ; i++ ) out_re[i] = in_re[i] * pow( -1.0, i ) ;
	}
	else
	{
		for( int i = 0 ; i < dimX ; i++ )
		{
			out_re[i] = in_re[i] * pow( -1.0, i ) ;
			out_im[i] = in_im[i] * pow( -1.0, i ) ;
		}
	}
}



void shift1d( int dimX, double * in_re, double * in_im, double * out_re, double * out_im )
{
	if( in_im == NULL )
	{
		for( int i = 0 ; i < dimX ; i++ ) out_re[i] = in_re[i] * pow( -1.0, i ) ;
	}
	else
	{
		for( int i = 0 ; i < dimX ; i++ )
		{
			out_re[i] = in_re[i] * pow( -1.0, i ) ;
			out_im[i] = in_im[i] * pow( -1.0, i ) ;
		}
	}
}



void shift2d( int dimX, int dimY, float  * in_re, float  * in_im, float  * out_re, float  * out_im )
{
	if( in_im == NULL )
	{
		for( int j = 0 ; j < dimY ; j++ )
			for( int i = 0 ; i < dimX ; i++ )
				out_re[ j*dimX+i ] = in_re[ j*dimX+i ] * pow( -1.0, i+j ) ;
	}
	else
	{
		for( int j = 0 ; j < dimY ; j++ )
			for( int i = 0 ; i < dimX ; i++ )
			{
				out_re[ j*dimX+i ] = in_re[ j*dimX+i ] * pow( -1.0, i+j ) ;
				out_im[ j*dimX+i ] = in_im[ j*dimX+i ] * pow( -1.0, i+j ) ;
			}
	}
}



void shift2d( int dimX, int dimY, double * in_re, double * in_im, double * out_re, double * out_im )
{
	if( in_im == NULL )
	{
		for( int j = 0 ; j < dimY ; j++ )
			for( int i = 0 ; i < dimX ; i++ )
				out_re[ j*dimX+i ] = in_re[ j*dimX+i ] * pow( -1.0, i+j ) ;
	}
	else
	{
		for( int j = 0 ; j < dimY ; j++ )
			for( int i = 0 ; i < dimX ; i++ )
			{
				out_re[ j*dimX+i ] = in_re[ j*dimX+i ] * pow( -1.0, i+j ) ;
				out_im[ j*dimX+i ] = in_im[ j*dimX+i ] * pow( -1.0, i+j ) ;
			}
	}
}



void shift3d( int dimX, int dimY, int dimZ, float  * in_re, float  * in_im, float  * out_re, float  * out_im )
{
	if( in_im == NULL )
	{
		for( int k = 0 ; k < dimZ ; k++ )
			for( int j = 0 ; j < dimY ; j++ )
				for( int i = 0 ; i < dimX ; i++ )
					out_re[ k*dimY*dimX+j*dimX+i ] = in_re[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
	}
	else
	{
		for( int k = 0 ; k < dimZ ; k++ )
			for( int j = 0 ; j < dimY ; j++ )
				for( int i = 0 ; i < dimX ; i++ )
				{
					out_re[ k*dimY*dimX+j*dimX+i ] = in_re[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
					out_im[ k*dimY*dimX+j*dimX+i ] = in_im[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
				}
	}
}



void shift3d( int dimX, int dimY, int dimZ, double * in_re, double * in_im, double * out_re, double * out_im )
{
	if( in_im == NULL )
	{
		for( int k = 0 ; k < dimZ ; k++ )
			for( int j = 0 ; j < dimY ; j++ )
				for( int i = 0 ; i < dimX ; i++ )
					out_re[ k*dimY*dimX+j*dimX+i ] = in_re[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
	}
	else
	{
		for( int k = 0 ; k < dimZ ; k++ )
			for( int j = 0 ; j < dimY ; j++ )
				for( int i = 0 ; i < dimX ; i++ )
				{
					out_re[ k*dimY*dimX+j*dimX+i ] = in_re[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
					out_im[ k*dimY*dimX+j*dimX+i ] = in_im[ k*dimY*dimX+j*dimX+i ] * pow( -1.0, i+j+k ) ;
				}
	}
}



void shift( int rank, int dimX, int dimY, int dimZ, float * in_re, float * in_im, float * out_re, float * out_im )
{
	switch( rank )
	{
		case 1:
			shift1d( dimX, in_re, in_im, out_re, out_im ) ;
			break ;
						
		case 2:
			shift2d( dimX, dimY, in_re, in_im, out_re, out_im ) ;
			break ;
						
		case 3:
			shift3d( dimX, dimY, dimZ, in_re, in_im, out_re, out_im ) ;
			break ;
					
		default:
			break ;
	}
}



void shift( int rank, int dimX, int dimY, int dimZ, double * in_re, double * in_im, double * out_re, double * out_im )
{
	switch( rank )
	{
		case 1:
			shift1d( dimX, in_re, in_im, out_re, out_im ) ;
			break ;
						
		case 2:
			shift2d( dimX, dimY, in_re, in_im, out_re, out_im ) ;
			break ;
						
		case 3:
			shift3d( dimX, dimY, dimZ, in_re, in_im, out_re, out_im ) ;
			break ;
					
		default:
			break ;
	}
}
