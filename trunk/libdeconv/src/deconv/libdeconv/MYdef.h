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
 * Filename:  MYdef.h
 */


#ifndef MYDEF_H
#define MYDEF_H


#include <boost/shared_array.hpp>
#include <cstdio>

typedef boost::shared_array< double >          DoubleArray ;
typedef boost::shared_array< float  >          SingleArray ;
typedef boost::shared_array< int    >          IntArray ;
typedef boost::shared_array< unsigned short >  ShortArray ;
typedef boost::shared_array< unsigned char >   ByteArray ;


#endif /* include MYdef.h */

