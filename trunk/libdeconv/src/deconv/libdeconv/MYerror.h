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
 * Filename:  MYerror.h
 */


#ifndef MYERROR_H
#define MYERROR_H


#include <errno.h>
#include <string.h>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/shared_ptr.hpp>


class Error : public std::exception
{
	public: 
	typedef boost::shared_ptr< Error > Ptr ;

	
	Error() : std::exception() {}


	Error( const Error& e ) : std::exception()
	{
		_error.str( e._error.str() ) ;
	}


	Error( std::string e ) : std::exception()
	{
		_error << e ;
	}
	

	virtual ~Error() throw() {}


	static Ptr create()
	{
		Ptr ret( new Error() ) ;
		return ret ;
	}


	static Ptr create( const Error& e )
	{
		Ptr ret( new Error( e ) ) ;
		return ret ;
	}		


	static Ptr create( std::string e )
	{
		Ptr ret( new Error( e ) ) ;
		return ret ;
	}		


	char const* what() const throw()
	{
		return _error.str().c_str() ;
	}

	
	Error& operator=( const Error& e )
	{
		_error.str( e._error.str() ) ;
		return *this ;
	}

	
	template< typename T >
	std::ostream& operator<<( T val )
	{
		return ( _error << val ) ;
	}


	std::ostream& operator<<( std::ostream& (*__pf) ( std::ostream& ) )
	{
		return ( _error << __pf ) ;
	}
	
       
	protected:
	std::stringstream  _error ;
} ;



class ErrnoError : public Error
{
	public:
	ErrnoError( std::string error )
	{
		int err = errno ;
		_error << strerror( err )
		       << " :: "
		       << error ;
	}
} ;



class ReadDataError : public Error
{
	public:
	ReadDataError( std::string s )
	{
		_error << " Error : failed to read data from " << s << ".\n" ;
	}
} ;



class WriteDataError : public Error
{
	public:
	WriteDataError( std::string s )
	{
		_error << " Error : failed to write data to " << s << ".\n" ;
	}
} ;


#endif /* include "MYerror.h" */
		
