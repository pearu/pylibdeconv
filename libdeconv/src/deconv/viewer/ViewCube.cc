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
 * Filename:  ViewCube.cc
 */


/*
 *	ViewCube is designed for viewing and handling a cubic data set.
 *
 *	One input argument is required by the program: <cube_name>
 *	For example : $./ViewCube MyCube.u8
 */


//#include <cmath>
#include <GL/glut.h>
#include "libdeconv/CCube.h"
#include "libdeconv/FFTW3fft.h"
#include "ViewerUtils.h"
#include "gnuplot_i.h"


#define WINDOW_WIDTH           512
#define WINDOW_HEIGHT          512
#define WINDOW_UPPER_LEFT_X    250
#define WINDOW_UPPER_LEFT_Y    250


static int       MainWindowID = 0 ;
static int       texImageX = 0 ;
static int       texImageY = 0 ;
static int       texImageZ = 0 ;
static unsigned char* texImageXY ;
static unsigned char* texImageXZ ;
static unsigned char* texImageYZ ;
static GLuint    texName[3] ;

static gnuplot_ctrl * PlotWindow ;
static bool           plot_refresh_on = true ;
static bool           plot_txtfile_on = false ;

static DoubleSlice  slice ;
static DoubleCube   cube ;
static DoubleCube   cubecopy ;
static char         cube_name[256] ;
static char         cube_head[256] ;
static int          cube_type = 0 ;
static std::string  cube_type_describe[4] = { "8-bit unsigned char", "16-bit unsigned short", 
                                              "32-bit single floating real", "64-bit double floating real" } ;

static double IntensityCenter = 0.5 ;
static double IntensityWindow = 1.0 ;

static float  spinX = 0 ;
static float  spinZ = 0 ;
static int    pageX = 0 ;
static int    pageY = 0 ;
static int    pageZ = 0 ;
static bool   XYviewON = true ;
static bool   XZviewON = true ;
static bool   YZviewON = true ;
static bool   XaxisON = true ;
static bool   YaxisON = true ;
static bool   ZaxisON = true ;
static bool   XcrossON = false ;
static bool   YcrossON = false ;
static bool   ZcrossON = false ;
static bool   labelON = true ;
static bool   XYlabelON = true ;
static bool   XZlabelON = true ;
static bool   YZlabelON = true ;
static float  labelX = 0.0 ;
static float  labelY = 0.0 ;
static float  labelZ = 0.0 ;

static bool   PowerSpectraON = false ;


std::string usage = "$./ViewCube <cube_name>\n" ;


void gnu_plot_data( double * data, int ndata, char * title )
{
	if( plot_txtfile_on )
	{
		std::string filename = std::string( title ) + ".txt" ;
		FILE * fp = fopen( filename.c_str(), "wt" ) ;
		if( fp )
		{
			for( int i = 0 ; i < ndata ; i++ ) fprintf( fp, "%d %f\n", i, data[i] ) ;
			fclose( fp ) ;
		}
		else throw ErrnoError( filename ) ;
	}
		
	double* index =  new double[ndata] ;
  	for ( int i = 0 ; i < ndata ; i++ ) index[i] = (double) i ;

	if( plot_refresh_on ) gnuplot_resetplot( PlotWindow ) ;
  	gnuplot_setstyle( PlotWindow, "lines" ) ;
  	gnuplot_plot_xy( PlotWindow, index, data, ndata, title ) ;
	delete [] index ;
}



void print_help( void )
{
	std::cout << "\n\n" ;

	std::cout << " ViewCube Menu (see instruction in README).\n" ;
	std::cout << " ViewCube Functional Keys :\n"  ;
	std::cout << " i, j, k, l  - rotate the object in the viewer\n" ;
	std::cout << " r           - reset viewer settings\n" ;
	std::cout << " q, a        - set sliceYZ    index    +1(q) and    -1(a)\n"
          	  << "             - move  label position  left(q) and right(a)\n" ;
	std::cout << " w, s        - set sliceXZ    index    +1(w) and    -1(s)\n"
        	  << "             - move  label position    up(w) and  down(s)\n" ;
  	std::cout << " e, d        - set sliceXY    index    +1(e) and    -1(d)\n"
  	          << "             - move  label position ahead(e) and  back(d)\n" ;
  	std::cout << " f           - swithch control functional keys : q,a,w,s,e,d\n" ;
 	std::cout << " p           - output the current viewer screen to a ppm file\n" ; 
  	std::cout << " h           - print help.\n" ;
  	std::cout << " ESC         - quit viewer\n" ;  	
  	std::cout << "\n\n" ;
}



void print_info( void )
{
  	std::cout << "\n\n" ;
  	std::cout << " Intensity_Center>>> " << IntensityCenter << "\n" ;
  	std::cout << " Intensity_Window>>> " << IntensityWindow << "\n" ;
  	std::cout << "        cube_name>>> " << cube_name << "\n" ;
  	std::cout << "        cube_type>>> " << cube_type_describe[cube_type-1] << "\n" ;
  	std::cout << "     PowerSpectra>>> " << (int) PowerSpectraON << "\n" ;
	if( XYviewON && XZviewON && YZviewON )
	{
		std::cout << "           Length>>> " << cube.length() << "\n" ;
		std::cout << "            Width>>> " << cube.width() << "\n" ;
  		std::cout << "           Height>>> " << cube.height() << "\n" ;
  		std::cout << "              Max>>> " << cube.cubemaxval() << "\n" ;
 		std::cout << "              Min>>> " << cube.cubeminval() << "\n" ;
  		std::cout << "             Mean>>> " << cube.cubemeanval() << "\n" ;
	}
	else
	{
		if( XYviewON )
		{
			cube.getslice( pageZ, XY_SLICE, slice ) ;
			std::cout << "          SliceXY>>> " << pageZ << "\n" ;		
		}
		if( XZviewON )
		{
			cube.getslice( pageY, XZ_SLICE, slice ) ;
			std::cout << "          SliceXZ>>> " << pageY << "\n" ;		
		}
		if( YZviewON )
		{
			cube.getslice( pageX, YZ_SLICE, slice ) ;
			std::cout << "          SliceYZ>>> " << pageX << "\n" ;		
		}		
		std::cout << "            Width>>> " << slice.width() << "\n" ;
  		std::cout << "           Height>>> " << slice.height() << "\n" ;
  		std::cout << "              Max>>> " << slice.slicemaxval() << "\n" ;
 		std::cout << "              Min>>> " << slice.sliceminval() << "\n" ;
  		std::cout << "             Mean>>> " << slice.slicemeanval() << "\n" ;
	}			
        std::cout << "\n\n" ;
}



void write_frame()
{
  	int       width  = glutGet( GLUT_WINDOW_WIDTH );
  	int       height = glutGet( GLUT_WINDOW_HEIGHT );  	
  	unsigned char* pixels = new unsigned char [ width * height * 3 ] ;

  	glReadBuffer( GL_BACK ) ;
  	glReadPixels( 0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels ) ;

  	write_my_byte_ppm( cube_head, width, height, pixels) ;
}



void write_slices( void )
{
	char setname[256], set_name[256], filehead[256] ;
	int  dir = 0, Dim1 = 0, Dim2 = 0, slice_plane = 0, EndSliceID = 0 ;
  	
	bool IsByte = true ;
	std::cout << " Indicate to save slices as 1-ByteImage or 0-ShortImage : " ;
	std::cin  >> dir ;
	IsByte = (bool) dir ;
  
	std::cout << " Indicate to save slices along 1-X or 2-Y or 3-Z : " ;
	std::cin  >> dir ;

	std::cout << " Saved slices along X will be named as <setname_yz###.pgm>.\n" ;
	std::cout << " Saved slices along Y will be named as <setname_xz###.pgm>.\n" ;
	std::cout << " Saved slices along Z will be named as <setname_xy###.pgm>.\n" ;
 	std::cout << " Input setname <<<< " ;
	std::cin  >> setname ;
	
  	switch( dir )
	{
		case 1 :
			EndSliceID  = cube.length() ;
			slice_plane = YZ_SLICE ;  
			Dim1 = cube.width() ;
			Dim2 = cube.height() ;
			sprintf( set_name, "%s_yz", setname ) ;
			break ;

		case 2 :
			EndSliceID  = cube.width() ;
			slice_plane = XZ_SLICE ;
			Dim1 = cube.length() ;
			Dim2 = cube.height() ;
			sprintf( set_name, "%s_xz", setname ) ;
			break ;

		case 3 :
			EndSliceID  = cube.height() ;
			slice_plane = XY_SLICE ;
			Dim1 = cube.length() ;
			Dim2 = cube.width() ;
			sprintf( set_name, "%s_xy", setname ) ;
			break ;
    	     	     
		default:
			std::cout << " Wrong Input = " << dir << " : slices must be saved along X(1)/Y(2)/Z(3).\n" ;
			break ;
	}
  	
	if( dir == 1 || dir == 2 || dir == 3 )
	{
		if( IsByte )
		{
			ByteSlice bslice( Dim1, Dim2 ) ;
			for( int i = 0 ; i < EndSliceID ; i++ ) 
			{
				cube.getslice( i, slice_plane, slice ) ;      				
				for( int j = 0 ; j < slice.size() ; j++ ) 
				{
					(bslice.data())[j] = byte_contrast( (slice.data())[j], IntensityWindow, IntensityCenter ) ;
				}
				name_seq_files( filehead, set_name, i+1, EndSliceID ) ;
				bslice.write( filehead ) ;
			}
		}
		else
		{
			ShortSlice sslice( Dim1, Dim2 ) ;
			for( int i = 0 ; i < EndSliceID ; i++ ) 
			{
				cube.getslice( i, slice_plane, slice ) ;
				for( int j = 0 ; j < slice.size() ; j++ ) 
				{
					(sslice.data())[j] = short_contrast( (slice.data())[j], IntensityWindow, IntensityCenter ) ;
				}
				name_seq_files( filehead, set_name, i+1, EndSliceID ) ;
				sslice.write( filehead ) ;
			}
		}	
	}
}



void write_cube( void )
{
	char filehead[256] ;
	int  saved_type ;
  	
	std::cout << " Input cube_name to be saved without suffix <<<< " ;
	std::cin  >> filehead ;
	for( int i = 0 ; i < 4 ; i++ ) 
	{
		std::cout << " cube_type = " << i+1 << " : " << cube_type_describe[i] << "\n" ;
	}
	std::cout << " Input cube_type to be saved <<<< " ;
	std::cin  >> saved_type ;
  	
	if( saved_type == 1 )
	{
		ByteCube bcube( cube.length(), cube.width(), cube.height() ) ;
		for( int i = 0 ; i < cube.size() ; i++ ) 
		{
			(bcube.data())[i] = byte_contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ;
		}  			
		bcube.write( filehead ) ;
	} 
	else if( saved_type == 2 )
	{
		ShortCube scube( cube.length(), cube.width(), cube.height() ) ; 
		for( int i = 0 ; i < cube.size() ; i++ ) 
		{
			(scube.data())[i] = short_contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ;
		}
		scube.write( filehead ) ;
	}
	else if( saved_type == 3 )
	{
		SingleCube fcube( cube.length(), cube.width(), cube.height() ) ;
		for( int i = 0 ; i < cube.size() ; i++ ) 
		{
			(fcube.data())[i] = (float) contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ;
 		}
		fcube.write( filehead ) ;
	}
	else if( saved_type == 4 )
	{
		DoubleCube dcube( cube.length(), cube.width(), cube.height() ) ;
		for( int i = 0 ; i < cube.size() ; i++ ) 
		{
			(dcube.data())[i] = contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ;
		}
		dcube.write( filehead ) ;
	}
	else
	{
		std::cout << " Wrong Input = " << saved_type << " : cube_type must be 1/2/3/4.\n" ;
	}
}



void destroy( void )
{
	gnuplot_close( PlotWindow ) ;
	glutDestroyWindow( MainWindowID ) ;
	glDeleteTextures( 3, texName ) ;
	exit( 0 ) ;
}



void initDisplay( void )
{	
	texImageX = cube.length() ;
	texImageY = cube.width() ;
	texImageZ = cube.height() ;
	while ( IsNotPowerOf2 ( texImageX ) ) texImageX++ ;
	while ( IsNotPowerOf2 ( texImageY ) ) texImageY++ ;
	while ( IsNotPowerOf2 ( texImageZ ) ) texImageZ++ ;
	texImageXY =  new unsigned char [ texImageX*texImageY*3 ] ;
	texImageXZ = new unsigned char [ texImageX*texImageZ*3 ]  ;
	texImageYZ = new unsigned char [ texImageY*texImageZ*3 ];
	
	IntensityCenter = ( cube.cubemaxval() + cube.cubeminval() ) / 2.0 ;
	IntensityWindow = cube.cubemaxval() - cube.cubeminval() ;
	if( IntensityWindow == 0.0 ) 
	{
		if( cube.cubemaxval() == 0.0 ) IntensityWindow = 1.0 ;
		else                           IntensityWindow = cube.cubemaxval() ;
		IntensityCenter = IntensityWindow / 2.0 ;
	}
	
	spinX     = 0.0 ;
	spinZ     = 0.0 ;  	
	pageX     = (int) ( cube.length() / 2 ) ;
	pageY     = (int) ( cube.width()  / 2 ) ;
	pageZ     = (int) ( cube.height() / 2 ) ;
	XYviewON  = true ;
	XZviewON  = true ;
	YZviewON  = true ;
	XaxisON   = true ;
	YaxisON   = true ;
	ZaxisON   = true ;
	XcrossON  = false ;
	YcrossON  = false ;
	ZcrossON  = false ;
	labelON   = false ;
	XYlabelON = false ;
	XZlabelON = false ;
	YZlabelON = false ;
	labelX    = 0.0 ;
	labelY    = 0.0 ;
	labelZ    = 0.0 ;
}



void createTexture( unsigned char* texImage, int tex_width, int tex_height, int shift_width, int shift_height, int num )
{
	int           k ;
	unsigned char texdata ;
	
	for( int i = 0 ; i < tex_width*tex_height*3 ; i++ ) texImage[i] = (unsigned char) 0 ;
	for( int j = 0 ; j < slice.height() ; j++ )
		for( int i = 0 ; i < slice.width() ; i++ ) 
		{
			texdata = byte_contrast( slice(i, j), IntensityWindow, IntensityCenter ) ;
			k = ( i+shift_width + (j+shift_height)*tex_width ) * 3 ;
			texImage[k]   = texdata ;
			texImage[k+1] = texdata ;
			texImage[k+2] = texdata ;
		}
  		
	glGenTextures( 1, &texName[num] ) ;
	glBindTexture( GL_TEXTURE_2D, texName[num] ) ;
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT ) ;
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT ) ;
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST ) ;
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST ) ;
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, tex_width, tex_height, 0, GL_RGB, GL_UNSIGNED_BYTE, texImage ) ;
}



void makeTexture()
{
	int texshiftx = (int) ( ( texImageX - cube.length() ) / 2 );
	int texshifty = (int) ( ( texImageY - cube.width()  ) / 2 );
	int texshiftz = (int) ( ( texImageZ - cube.height() ) / 2 );
	cube.getslice( pageZ, XY_SLICE, slice ) ;
	createTexture ( texImageXY, texImageX, texImageY, texshiftx, texshifty, 0 ) ;	
	cube.getslice( pageY, XZ_SLICE, slice ) ;
	createTexture ( texImageXZ, texImageX, texImageZ, texshiftx, texshiftz, 1 ) ;
	cube.getslice( pageX, YZ_SLICE, slice ) ;
	createTexture ( texImageYZ, texImageY, texImageZ, texshifty, texshiftz, 2 ) ;
}



void display( void )
{
	char  label[256] ;
	float vx = (float) cube.length() ;
	float vy = (float) cube.width() ;
	float vz = (float) cube.height() ;
	float texvx = vx / ( (float) texImageX ) ;
	float texvy = vy / ( (float) texImageY ) ;
	float texvz = vz / ( (float) texImageZ ) ;
	float vmax = vx ;
	if ( vmax < vy  ) vmax = vy ;
	if ( vmax < vz  ) vmax = vz ;
	vx = vx / vmax ;
	vy = vy / vmax ;
	vz = vz / vmax ;

	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glShadeModel( GL_FLAT ) ;
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT) ;
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1) ;
	glPushMatrix() ;
	glRotatef( spinZ, 0.0, 1.0, 0.0 ) ;
	glRotatef( spinX, 1.0, 0.0, 0.0 ) ;
	glEnable( GL_DEPTH_TEST ) ;
	glEnable( GL_TEXTURE_2D ) ;
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE ) ;
  	
	if( XYviewON ) 
	{
		glBindTexture( GL_TEXTURE_2D, texName[0] ) ;
		glBegin( GL_QUADS ) ;
		  glTexCoord2f( (1.0-texvx)/2.0, (1.0-texvy)/2.0 ) ;  glVertex3f( -vx, -vz, -vy) ;    // Back  Left  (Bottom)
		  glTexCoord2f( (1.0-texvx)/2.0, (1.0+texvy)/2.0 ) ;  glVertex3f( -vx, -vz,  vy) ;    // Front Left  (Bottom)
		  glTexCoord2f( (1.0+texvx)/2.0, (1.0+texvy)/2.0 ) ;  glVertex3f(  vx, -vz,  vy) ;    // Front Right (Bottom)
		  glTexCoord2f( (1.0+texvx)/2.0, (1.0-texvy)/2.0 ) ;  glVertex3f(  vx, -vz, -vy) ;    // Back  Right (Bottom)
		glEnd() ;
	}
  	
	if( XZviewON ) 
	{
		glBindTexture( GL_TEXTURE_2D, texName[1] ) ;
		glBegin(GL_QUADS) ;
		  glTexCoord2f( (1.0-texvx)/2.0, (1.0-texvz)/2.0 ) ;  glVertex3f( -vx, -vz, -vy ) ;    // Top    Left  (Back)
		  glTexCoord2f( (1.0-texvx)/2.0, (1.0+texvz)/2.0 ) ;  glVertex3f( -vx,  vz, -vy ) ;    // Bottom Left  (Back)
		  glTexCoord2f( (1.0+texvx)/2.0, (1.0+texvz)/2.0 ) ;  glVertex3f(  vx,  vz, -vy ) ;    // Bottom Right (Back)
		  glTexCoord2f( (1.0+texvx)/2.0, (1.0-texvz)/2.0 ) ;  glVertex3f(  vx, -vz, -vy ) ;    // Top    Right (Back)
		glEnd();
	}
  	
	if( YZviewON ) 
	{
		glBindTexture( GL_TEXTURE_2D, texName[2] ) ;
		glBegin(GL_QUADS) ;
		  glTexCoord2f( (1.0-texvy)/2.0, (1.0-texvz)/2.0 ) ;  glVertex3f( -vx, -vz, -vy ) ;    // Front Top    (Left)
		  glTexCoord2f( (1.0-texvy)/2.0, (1.0+texvz)/2.0 ) ;  glVertex3f( -vx,  vz, -vy ) ;    // Front Bottom (Left)
		  glTexCoord2f( (1.0+texvy)/2.0, (1.0+texvz)/2.0 ) ;  glVertex3f( -vx,  vz,  vy ) ;    // Back  Bottom (Left)
		  glTexCoord2f( (1.0+texvy)/2.0, (1.0-texvz)/2.0 ) ;  glVertex3f( -vx, -vz,  vy ) ;    // Back  Top    (Left)
		glEnd();
	}
  	
	glDisable( GL_TEXTURE_2D ) ;
	glDisable( GL_DEPTH_TEST ) ;

	if( XaxisON ) 
	{
		glColor3f( 1.0, 0.0, 0.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, -vz, -vy ) ;
		  glVertex3f(  vx, -vz, -vy ) ;
		glEnd() ;
		glRasterPos3f( vx*1.1, -vz, -vy ) ;
		glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, 'X' ) ;
	}

	if( YaxisON ) 
	{
		glColor3f( 0.0, 1.0, 0.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, -vz, -vy ) ;
		  glVertex3f( -vx, -vz,  vy ) ;
		glEnd() ;
		glRasterPos3f( -vx, -vz, vy*1.1 ) ;
		glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, 'Y' ) ;
	}

	if( ZaxisON ) 
	{
		glColor3f( 0.0, 0.0, 1.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, -vz, -vy ) ;
		  glVertex3f( -vx,  vz, -vy ) ;
		glEnd() ;
		glRasterPos3f( -vx, 1.1*vz, -vy ) ;
		glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, 'Z' ) ;
	}

	if( XcrossON ) 
	{
		glColor3f( 0.0, 1.0, 0.0 ) ;
		glBegin(GL_LINES) ;
		  glVertex3f( ((float)(pageX-cube.length()/2))/((float)(cube.length()/2))*vx, -vz, -vy ) ;
		  glVertex3f( ((float)(pageX-cube.length()/2))/((float)(cube.length()/2))*vx, -vz,  vy ) ;
		glEnd() ;
		glColor3f( 0.0, 0.0, 1.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( ((float)(pageX-cube.length()/2))/((float)(cube.length()/2))*vx, -vz, -vy ) ;
		  glVertex3f( ((float)(pageX-cube.length()/2))/((float)(cube.length()/2))*vx,  vz, -vy ) ;
		glEnd() ;
	}

	if( YcrossON ) 
	{
		glColor3f( 1.0, 0.0, 0.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, -vz, ((float)(pageY-cube.width()/2))/((float)(cube.width()/2))*vy ) ;
		  glVertex3f(  vx, -vz, ((float)(pageY-cube.width()/2))/((float)(cube.width()/2))*vy ) ;
		glEnd() ;
		glColor3f( 0.0, 0.0, 1.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, -vz, ((float)(pageY-cube.width()/2))/((float)(cube.width()/2))*vy ) ;
		  glVertex3f( -vx,  vz, ((float)(pageY-cube.width()/2))/((float)(cube.width()/2))*vy ) ;
		glEnd() ;
	}

	if( ZcrossON ) 
	{
		glColor3f( 1.0, 0.0, 0.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, ((float)(pageZ-cube.height()/2))/((float)(cube.height()/2))*vz, -vy ) ;
		  glVertex3f(  vx, ((float)(pageZ-cube.height()/2))/((float)(cube.height()/2))*vz, -vy ) ;
		glEnd() ;
		glColor3f( 0.0, 1.0, 0.0 ) ;
		glBegin( GL_LINES ) ;
		  glVertex3f( -vx, ((float)(pageZ-cube.height()/2))/((float)(cube.height()/2))*vz, -vy ) ;
		  glVertex3f( -vx, ((float)(pageZ-cube.height()/2))/((float)(cube.height()/2))*vz,  vy ) ;
		glEnd();
	}

	glPopMatrix() ;

	if( XYlabelON || XZlabelON || YZlabelON ) 
	{
		glColor3f( 1.0, 1.0, 1.0 ) ;    		
		if ( XYlabelON && XZlabelON && YZlabelON )
		{
			sprintf( label, "SliceXY = %d/%d  SliceXZ = %d/%d  SliceYZ = %d/%d",
			pageZ+1, cube.height(), pageY+1, cube.width(), pageX+1, cube.length() ) ;
		}
		else if( XYlabelON && XZlabelON && (!YZlabelON) )
		{
			sprintf( label, "SliceXY = %d/%d  SliceXZ = %d/%d", pageZ+1, cube.height(), pageY+1, cube.width() ) ;
		}   		
		else if( XYlabelON && YZlabelON && (!XZlabelON) )
		{
			sprintf( label, "SliceXY = %d/%d  SliceYZ = %d/%d", pageZ+1, cube.height(), pageX+1, cube.length() ) ;
		}
		else if( XZlabelON && YZlabelON && (!XYlabelON) )
		{
			sprintf( label, "SliceXZ = %d/%d  SliceYZ = %d/%d", pageY+1, cube.width(), pageX+1, cube.length() ) ;
		}
		else if( YZlabelON && (!XZlabelON) && (!XYlabelON) )
		{ 
			sprintf( label, "SliceYZ = %d/%d", pageX+1, cube.length() ) ;
		}
		else if( XZlabelON && (!YZlabelON) && (!XYlabelON) )
		{ 
			sprintf( label, "SliceXZ = %d/%d", pageY+1, cube.width() ) ;
		}
		else    
		{
			sprintf( label, "SliceXY = %d/%d", pageZ+1, cube.height() ) ;
		}
		glRasterPos3f( labelX, labelY, labelZ ) ;
		for( int i = 0 ; i < ((int) strlen( label )) ; i++ ) glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, label[i] ) ; 
	}

	glFlush() ;
	glutSwapBuffers() ;
}



void reshape( int w, int h )
{
	glViewport( 0, 0, (GLsizei)w, (GLsizei)h ) ;
	glMatrixMode( GL_PROJECTION ) ;
	glLoadIdentity() ;
	gluPerspective( 60.0, (GLfloat) w/(GLfloat) h, 1.0, 30.0 ) ;
	glMatrixMode( GL_MODELVIEW ) ;
	glLoadIdentity() ;
	glTranslatef( 0.0, 0.0, -3.6 ) ;  	
}



void keyboard( unsigned char key, int x, int y )
{
	switch( key ) 
	{
		case 108: // l - rotate object
			spinZ = spinZ + 2.0 ; 
			glutPostRedisplay() ;
    		     	break ;

		case 106: // j - rotate object 
			spinZ = spinZ - 2.0 ;
			glutPostRedisplay() ;
			break ;

		case 105: // i - rotate object
			spinX = spinX + 2.0 ;
			glutPostRedisplay() ;
			break ;

		case 107: // k - rotate object
    		     	spinX = spinX - 2.0 ;
    		     	glutPostRedisplay() ;
    		     	break ;

		case 114: // r - reset  viewer
			initDisplay() ;
			makeTexture() ;
			glutPostRedisplay() ;
			break ;

		case 113: // q - pageX up or labelX left 
			if( labelON )
			{
				labelX = labelX - 0.1 ;
			}
			else 
			{
				if( pageX == cube.length()-1 ) pageX = 0 ;
				else                           pageX++ ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 97: // a - pageX down or labelX right
			if( labelON ) 
			{
				labelX = labelX + 0.1 ;
			}
			else 
			{
				if( pageX == 0 ) pageX = cube.length() - 1 ;
				else             pageX-- ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 119: // w - pageY up or labelY up
			if( labelON ) 
			{
				labelY = labelY + 0.1 ;
			}
			else 
			{
				if( pageY == cube.width()-1 ) pageY = 0 ;
				else                          pageY++ ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 115: // s - pageY down or labelY down
			if( labelON ) 
			{
				labelY = labelY - 0.1 ;
			}
			else 
			{
				if( pageY == 0 ) pageY = cube.width() - 1 ;
				else             pageY-- ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 101: // e - pageZ up or labelZ forward 
			if( labelON )
			{ 
				labelZ = labelZ - 0.1 ;
			}
			else 
			{
				if( pageZ == cube.height()-1 ) pageZ = 0 ;
				else                           pageZ++ ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 100: // d - pageZ down or labelZ backward
			if( labelON ) 
			{
				labelZ = labelZ + 0.1 ;
			}
			else 
			{
				if( pageZ == 0 ) pageZ = cube.height() - 1 ;
				else             pageZ-- ;
				makeTexture() ;
			}
			glutPostRedisplay() ;
			break ;

		case 102: // f - switch functional keys
			if( labelON ) labelON = false ;
			else          labelON = ( XYlabelON || XZlabelON || YZlabelON ) ;
			break ;

		case 112: // p - snap current viewer-screen frame
			write_frame() ;
			break ;

		case 104: // h - help
			print_help() ;
			break ;

  		case 27: // ESC - quit viewer
			destroy();
			break ;

		default:
			break ;
	}
}



void mouse( int button, int state, int x, int y )
{
	GLubyte curpixel[3] ;
	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) 
	{
		glReadPixels( x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, curpixel ) ;
		printf(" Screen Position: x=%4d , y=%4d -> Value(RGB): %3d , %3d , %3d\n",
		x ,y, curpixel[0], curpixel[1], curpixel[2]) ;
	}
}



void fileMenuFunc( int choice )
{
	switch( choice ) 
	{
		case 1:
			cube.read( cube_name ) ; 
			initDisplay() ;
			makeTexture() ;
			glutPostRedisplay() ;
			break ;

		case 2: 
			write_cube() ;
			break ;

		case 3: 
			write_slices() ;
			break;

		default:
			break;
	}
}



void viewMenuFunc( int choice )
{
	switch( choice ) 
	{
		case 1: // view XY-projection
			spinX    = -90 ;
			spinZ    = 0 ;
			XaxisON  = true ;
			YaxisON  = true ;
			ZaxisON  = false ;
			XYviewON = true ;
			XZviewON = false ;
			YZviewON = false ;
			break ;

		case 2: // view XZ-projection
			spinX    = 0 ;
			spinZ    = 180 ;
			XaxisON  = true ;
			YaxisON  = false ;
			ZaxisON  = true ;
			XYviewON = false ;
			XZviewON = true ;
			YZviewON = false ;
			break ;

		case 3: // view YZ-projection
			spinX    = 0 ;
			spinZ    = 90 ;
			XaxisON  = false ;
			YaxisON  = true ;
			ZaxisON  = true ;
			XYviewON = false ;
			XZviewON = false ;
			YZviewON = true ;
    	     		break ;

		case 4: // view all projections
			spinX    = 0 ;
			spinZ    = 0 ;
			XaxisON  = true ;
			YaxisON  = true ;
			ZaxisON  = true ;
			XYviewON = true ;
			XZviewON = true ;
			YZviewON = true ;
			break ;

		default:
			break;
	}
	glutPostRedisplay() ;
}



void plotMenuFunc( int choice )
{
	int         ndata = 0, slice_plane = 0, maxindex = 0 ;
	int         remove_peak = 0 ;
	double      maxdata = 0.0, sumdata = 0.0 ; 
	char        title[256] ;
	double*     data = NULL ;
	std::string s ;

	switch( choice ) 
	{
		case 1: 
		{
			ndata = 256 ;			
			data = new double[ ndata ]  ;
			for( int i = 0 ; i < ndata ; i++ ) data[i] = 0.0 ;
			if( XYviewON && XZviewON && YZviewON )
			{
				for( int i = 0 ; i < cube.size() ; i++ ) 
				{
					data[ (int) byte_contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ]++ ;
				}
				sprintf( title, "ByteHistogram_%s", cube_head ) ; 
			}
			else
			{	
				if( XYviewON ) 
				{
					cube.getslice( pageZ, XY_SLICE, slice ) ;
					sprintf( title, "ByteHistogram_%s_sliceXY%d", cube_head, pageZ ) ;
				}
				if( XZviewON )
				{
					cube.getslice( pageY, XZ_SLICE, slice ) ;
					sprintf( title, "ByteHistogram_%s_sliceXZ%d", cube_head, pageY ) ;
				}
				if( YZviewON )
				{
					cube.getslice( pageX, YZ_SLICE, slice ) ;
					sprintf( title, "ByteHistogram_%s_sliceYZ%d", cube_head, pageX );
				}
				for( int i = 0 ; i < slice.size() ; i++ )
				{
					data[ (int) byte_contrast( (slice.data())[i], IntensityWindow, IntensityCenter ) ]++ ;
				}
			}         				
			for( int i = 0 ; i < ndata ; i++ )
			{
				sumdata += data[i] ;
				if( data[i] >= maxdata ) 
				{ 
					maxdata = data[i] ; 
					maxindex = i ; 
				}
			}
			std::cout << " ByteHistogram peak found at " << maxindex << " : " << maxdata << "/" << sumdata << "\n" ;
			std::cout << " Indicate whether removing(1) the peak or not(0) to plot data <<<< " ;
			std::cin  >> remove_peak ;
			if( remove_peak ) data[ maxindex ] = 0.0 ;
			gnu_plot_data( data, ndata, title ) ;
			break ;
		}
    			
		case 2: 
		{
			ndata = 65536 ; 			
			data = new double[ ndata ] ;
			for( int i = 0 ; i < ndata ; i++ ) data[i] = 0.0 ;
			if( XYviewON && XZviewON && YZviewON )
			{
				for( int i = 0 ; i < cube.size() ; i++ ) 
				{
					data[ (int) short_contrast( (cube.data())[i], IntensityWindow, IntensityCenter ) ]++ ;
				}
				sprintf( title, "ShortHistogram_%s", cube_head ) ; 
			}
			else
			{	
				if( XYviewON ) 
				{
					cube.getslice( pageZ, XY_SLICE, slice ) ;
					sprintf( title, "ShortHistogram_%s_sliceXY%d", cube_head, pageZ ) ;
				}
				if( XZviewON ) 
				{
					cube.getslice( pageY, XZ_SLICE, slice ) ;
					sprintf( title, "ShortHistogram_%s_sliceXZ%d", cube_head, pageY ) ;
				}
				if( YZviewON ) 
				{
					cube.getslice( pageX, YZ_SLICE, slice ) ;                  			
					sprintf( title, "ShortHistogram_%s_sliceYZ%d", cube_head, pageX );
				}
				for( int i = 0 ; i < slice.size() ; i++ )
				{
					data[ (int) short_contrast( (slice.data())[i], IntensityWindow, IntensityCenter ) ]++ ;
				}
			}
			for( int i = 0 ; i < ndata ; i++ )
			{
				sumdata += data[i] ;
				if( data[i] >= maxdata ) 
				{ 
					maxdata = data[i] ; 
					maxindex = i ; 
				}
			}
			std::cout << " ShortHistogram peak found at " << maxindex << " : " << maxdata << "/" << sumdata << "\n" ;
			std::cout << " Indicate whether removing(1) the peak or not(0) to plot data <<<< " ;
			std::cin  >> remove_peak ;
			if( remove_peak ) data[ maxindex ] = 0.0 ;
			gnu_plot_data( data, ndata, title ) ;
			break ;
		}

		case 3:  
		{			 
			if( XYviewON ) 
			{
				ndata = cube.height() ;   
				slice_plane = XY_SLICE ;
				s = std::string( " Focus along Z found at SliceXY " ) ;
				sprintf( title, "FocusCurve_%s_alongZ", cube_head ) ;
			}    			
			if( XZviewON && (!YZviewON) && (!XYviewON) ) 
			{
				ndata = cube.width() ;
				slice_plane = XZ_SLICE ;
				s = std::string( " Focus along Y found at SliceXZ " ) ;
				sprintf( title, "FocusCurve_%s_alongY", cube_head ) ;
			}        		
			if( YZviewON && (!XZviewON) && (!XYviewON) ) 
			{
				ndata = cube.length() ;
				slice_plane = YZ_SLICE ;
				s = std::string( " Focus along X found at SliceYZ " ) ;
				sprintf( title, "FocusCurve_%s_alongX", cube_head ) ;
			}  
			data = new double[ ndata ] ;
			for( int k = 0 ; k < ndata ; k++ ) 
			{
				data[k] = 0.0 ;
				cube.getslice( k, slice_plane, slice ) ;
				for( int j = 0 ; j < slice.height() ; j++ )
					for( int i = 0 ; i < slice.width() ; i++ )
					{
						if( i < slice.width()-1 ) 
						{
							data[k] += ( slice(i,j) * ( slice(i,j)-slice(i+1,j) ) ) ;
						}
					}
				if( data[k] >= maxdata ) 
				{ 
					maxdata = data[k] ; 
					maxindex = k ; 
				}
			}
			std::cout << s << maxindex << " : " << maxdata << "\n" ;
			gnu_plot_data( data, ndata, title ) ;
			break;
		}
    			
		case 4:
		{
			if( XYviewON )
			{
				ndata = cube.height() ;
				slice_plane = XY_SLICE ;
				s = std::string( " Maximum slice-mean along Z found at SliceXY " ) ;
				sprintf( title, "MeanCurve_%s_alongZ", cube_head ) ;
			}			
			if( XZviewON && (!YZviewON) && (!XYviewON) ) 
			{
				ndata = cube.width() ;
				slice_plane = XZ_SLICE ;
				s = std::string( " Maximum slice-mean along Y found at SliceXZ " ) ;
				sprintf( title, "MeanCurve_%s_alongY", cube_head ) ;
			}	     	     			
			if( YZviewON && (!XYviewON) && (!XZviewON) )
			{
				ndata = cube.length() ;
				slice_plane = YZ_SLICE ;
				s = std::string( " Maximum slice-mean along X found at SliceYZ " ) ;  
				sprintf( title, "MeanCurve_%s_alongX", cube_head ) ;
			}   
			data = new double[ ndata ] ;   
			for( int i = 0 ; i < ndata ; i++ ) 
			{
				cube.getslice( i, slice_plane, slice ) ;
				data[i] = slice.slicemeanval() ; 
				if( data[i] >= maxdata ) 
				{ 
					maxdata = data[i] ; 
					maxindex = i ; 
				}
			}
			std::cout << s << maxindex << " : " << maxdata << "\n" ;
			gnu_plot_data( data, ndata, title ) ;
			break;
		}
    		
		case 5:
		{
			if( XYviewON )
			{
				ndata = cube.height() ;
				slice_plane = XY_SLICE ;
				s = std::string( " Maximum slice-max along Z found at SliceXY " ) ;
				sprintf( title, "MaxCurve_%s_alongZ", cube_head ) ;
			}			
			if( XZviewON && (!YZviewON) && (!XYviewON) ) 
			{
				ndata = cube.width() ;
				slice_plane = XZ_SLICE ;
				s = std::string( " Maximum slice-max along Y found at SliceXZ " ) ;
				sprintf( title, "MaxCurve_%s_alongY", cube_head ) ;
			}	     	     			
			if( YZviewON && (!XYviewON) && (!XZviewON) )
			{
				ndata = cube.length() ;
				slice_plane = YZ_SLICE ;
				s = std::string( " Maximum slice-max along X found at SliceYZ " ) ;  
				sprintf( title, "MaxCurve_%s_alongX", cube_head ) ;
			}
			data = new double[ ndata ]  ;	  
			for( int i = 0 ; i < ndata ; i++ ) 
			{
				cube.getslice( i, slice_plane, slice ) ;
				data[i] = slice.slicemaxval() ; 
				if( data[i] >= maxdata ) 
				{ 
					maxdata = data[i] ; 
 					maxindex = i ; 
				}
			}
			std::cout << s << maxindex << " : " << maxdata << "\n" ;
			gnu_plot_data( data, ndata, title ) ;
			break;
		}
    		
		case 6:
			plot_refresh_on = true ;
			break ;
    			
		case 7:
			plot_refresh_on = false ;
			break ;
    			
		case 8:
			plot_txtfile_on = true ;
			break ;
    		
		case 9:
			plot_txtfile_on = false ;
			break ;

		default:
			break ;
	}
	if (data)
		delete [] data ;
}



void operateMenuFunc( int choice )
{
	int    low, high, rx, ry, hz ;
	double value = 0.0 ;
	
	switch( choice ) 
	{
		case 1:
			std::cout << " Apply contrasting, input the new Intensity_Window_Width  <<<< " ;
			std::cin  >> value ;
			if( value <= 0.0 )
			{
				std::cout << " Wrong Input : Intensity_Window_Width must be positive.\n" ;
			}
			else
			{ 
				IntensityWindow = value ;
				std::cout << " Apply contrating, input the new Intensity_Window_Center <<<< " ;
				std::cin  >> IntensityCenter ;
				makeTexture() ;
				glutPostRedisplay() ;
			}
			break ;
    					
		case 2:
			std::cout << " Cube length (X) is indexed from 0 to " << cube.length()-1 << "\n" ;
			std::cout << " Cube width  (Y) is indexed from 0 to " << cube.width()-1  << "\n" ;
			std::cout << " Cube height (Z) is indexed from 0 to " << cube.height()-1 << "\n" ;
			std::cout << " Padding cube, input the padding value <<<< " ;
			std::cin  >> value ; 
			std::cout << " Padding cube to X-, input how many slices to be padded before X=0 <<<< " ;
			std::cin  >> low ;
			std::cout << " Padding cube to X+, input how many slices to be padded after X=" << cube.length()-1 << " <<<< " ;
			std::cin  >> high ;
			if( !cube.length_pad( low, high, value ) )
			{
				std::cout << " Padding cube to Y-, input how many slices to be padded before Y=0 <<<< " ;
				std::cin  >> low ;
				std::cout << " Padding cube to Y+, input how many slices to be padded after Y=" 
				          << cube.width()-1 << " <<<< " ;
				std::cin  >> high ;
				if( !cube.width_pad( low, high, value ) )
				{
					std::cout << " Padding cube to Z-, input how many slices to be padded before Z=0 <<<< " ;
					std::cin  >> low ;
					std::cout << " Padding cube to Z+, input how many slices to be padded after Z=" 
					          << cube.height()-1 << " <<<< " ;
					std::cin  >> high ;
					if( !cube.height_pad( low, high, value ) )
					{
						initDisplay() ;
						makeTexture() ;
						glutPostRedisplay() ;
					}
				}
			}
			break ;

		case 3:
			std::cout << " Cube length (X) is indexed from 0 to " << cube.length()-1 << "\n" ;
			std::cout << " Cube width  (Y) is indexed from 0 to " << cube.width()-1  << "\n" ;
			std::cout << " Cube height (Z) is indexed from 0 to " << cube.height()-1 << "\n" ;
			std::cout << " Cropping cube from X-, input the cropped lower  index : X <<<< " ;
			std::cin  >> low ;
			std::cout << " Cropping cube from X+, input the cropped higher index : X <<<< " ;
			std::cin  >> high ;
			if( !cube.length_crop( low, high ) )
			{
				std::cout << " Cropping cube from Y-, input the cropped lower  index : Y <<<< " ;
				std::cin  >> low ;
				std::cout << " Cropping cube from Y+, input the cropped higher index : Y <<<< " ;
				std::cin  >> high ;
				if( !cube.width_crop( low, high ) )
				{
					std::cout << " Cropping cube from Z-, input the cropped lower  index : Z <<<< " ;
					std::cin  >> low ;
 					std::cout << " Cropping cube from Z+, input the cropped higher index : Z <<<< " ;
					std::cin  >> high ;
					if( !cube.height_crop( low, high ) )
					{
						initDisplay() ;
						makeTexture() ;
						glutPostRedisplay() ;
					}
				}
			}
			break ;

		case 4:
			std::cout << " Apply substration : data_intensity = (data_intensity-value, 0)?(data_intensity-value>0).\n"
			          << " Now, please input the substracted value <<<< " ;
			std::cin  >> value ;	
			for( int i = 0 ; i < cube.size() ; i++ ) 
			{
				if( (cube.data())[i] - value < 0.0 ) (cube.data())[i] = 0.0 ;
				else                                 (cube.data())[i] = (cube.data())[i] - value ;
			}
			initDisplay() ;
			makeTexture() ;
			glutPostRedisplay() ;
			break ;
    			
		case 5: 
			if( !cube.shift() )
			{
				initDisplay() ;
				makeTexture() ;
				glutPostRedisplay() ;
			}
			break ;
    			
		case 6:
			std::cout << " Apply lnTransform : data_intensity = ln( data_intensity * 10^Power +1 ).\n" 
			          << " Now, please input the lnTransformed power <<<< " ;
			std::cin  >> low ;
			if( low >= 0 )
			{
				for( int i = 0 ; i < cube.size() ; i++ )
				{
					(cube.data())[i] = log( (cube.data())[i] * pow( 10.0, (double) low ) + 1.0 ) ;
				}
			}
			initDisplay() ;
			makeTexture() ;
			glutPostRedisplay() ;
			break ;

		case 7:
		{
			if( PowerSpectraON )
			{
				PowerSpectraON = false ;
				cube = cubecopy ;
				cubecopy.free() ;
			}
			else
			{
				PowerSpectraON = true ;
				cubecopy = cube ;
				double* buf = new double[cube.size()] ;
				for( int i = 0 ; i < cube.size() ; i++ ) buf[i] = 0.0 ;
				fft3d( cube.length(), cube.width(), cube.height(), true, true, 
				       cube.data(), buf, cube.data(), buf ) ;
				for( int i = 0 ; i < cube.size() ; i++ ) 
				{
					(cube.data())[i] = sqrt( (cube.data())[i] * (cube.data())[i] + buf[i] * buf[i] ) ;
				}
				delete [] buf ;
			}
			initDisplay() ;
			makeTexture() ;
			glutPostRedisplay() ;
			break ;
    		}

		case 8:
                {
                        ByteCube filter_cube( cube.length(), cube.width(), cube.height() ) ;
                        filter_cube.fillcube( (unsigned char) 0 ) ;
                        std::cout << " Apply a cylinder-shaped LowPass filter : \n" ;
                        std::cout << " Now, input cylinder radius along x-dimension <<<< " ;
                        std::cin  >> rx ;
                        std::cout << " Now, input cylinder radius along y-dimension <<<< " ;
                        std::cin  >> ry ;
                        std::cout << " Now, input cylinder height along z-dimension <<<< " ;
                        std::cin  >> hz ;
                        if( !filter_cube.draw_cylinder( cube.length()/2, cube.width()/2, cube.height()/2,
                                                        rx, ry, hz, (unsigned char) 1 ) )
                        {
								double* buf = new double[cube.size()] ;
                                for( int i = 0 ; i < cube.size() ; i++ ) buf[i] = 0.0 ;
                                fft3d( cube.length(), cube.width(), cube.height(), true, true,
                                       cube.data(), buf, cube.data(), buf ) ;
                                for( int i = 0 ; i < cube.size() ; i++ )
                                {
                                        (cube.data())[i] *= (double) (filter_cube.data())[i] ;
                                        buf[i] *= (double) (filter_cube.data())[i] ;
                                }
                                fft3d( cube.length(), cube.width(), cube.height(), false, true,
                                       cube.data(), buf, cube.data(), buf ) ;
                                for( int i = 0 ; i < cube.size(); i++ )
                                	if( (cube.data())[i] < 0 ) 
										(cube.data())[i] = 0.0 ;

								delete [] buf ;
                        }
                        initDisplay() ;
                        makeTexture() ;
                        glutPostRedisplay() ;
                }
		
		default:
			break ;
	}
}



void controlMenuFunc( int choice )
{
	switch( choice ) 
	{
		case 1: // x - axis
			if( XaxisON ) XaxisON = false ;
			else          XaxisON = true ;
			break ;

		case 2: // y - axis
			if( YaxisON ) YaxisON = false ;
			else          YaxisON = true ;
			break ;

		case 3: // z - axis
			if( ZaxisON ) ZaxisON = false ;
			else          ZaxisON = true ;
			break ;

		case 4: // x - cross
			if( XcrossON ) XcrossON = false ;
			else           XcrossON = true ;
			break ;

		case 5:	// y - cross
			if( YcrossON ) YcrossON = false ;
			else           YcrossON = true ;
			break ;

		case 6:	// z - cross
			if( ZcrossON ) ZcrossON = false ;
			else           ZcrossON = true ;
			break ;

		case 7:	// sliceXY - label
			if( XYlabelON ) XYlabelON = false ;
			else            XYlabelON = true ;
			break ;

		case 8:	// sliceXZ - label
			if( XZlabelON ) XZlabelON = false ;
			else            XZlabelON = true ;
			break ;

		case 9:	// sliceYZ - label
			if( YZlabelON ) YZlabelON = false ;
			else            YZlabelON = true ;
			break ;	
    	     		
		default:
			break ;
	}
	glutPostRedisplay() ;
}



void mainMenuFunc( int choice ) 
{
	switch( choice ) 
	{		
		case 1:
			print_help() ;
			break ;
	
		case 2:
			print_info() ;
			break ;

		case 3:
			destroy() ;
			break;

		default:
			break;
	}
}



int main(int argc, char** argv)
{
	if( argc != 2 ) 
	{ 
		std::cout << usage ;
		return 1 ; 
	} 
	else 
	{
		sprintf( cube_name, "%s", argv[1] ) ;
		cube_type = cube.read( cube_name ) ;
		std::string s = std::string( cube_name ) ;
		std::string::size_type pos0 = s.find_last_of( '/' ) ;
		std::string::size_type pos1 = s.find_last_of( '.' ) ;
		std::string::size_type pos_begin = pos0 + 1 ;
		std::string::size_type pos_end = pos1 - pos0 - 1 ;		
		if( pos0 == std::string::npos ) 
		{
			pos_begin = 0 ; 
			pos_end = pos1 ;
		}
		std::string ss = s.substr( pos_begin, pos_end ) ; 
		sprintf( cube_head, "%s", ss.c_str() ) ;

		PlotWindow = gnuplot_init() ;

		glutInit( &argc, argv ) ;
		glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB ) ;
		glutInitWindowSize( WINDOW_WIDTH, WINDOW_HEIGHT ) ;
		glutInitWindowPosition( WINDOW_UPPER_LEFT_X, WINDOW_UPPER_LEFT_Y ) ;
		MainWindowID = glutCreateWindow( cube_name ) ;

		print_help() ;
		initDisplay() ;
		makeTexture() ;

		glutDisplayFunc( display ) ;
		glutReshapeFunc( reshape ) ;
		glutKeyboardFunc( keyboard );
		glutMouseFunc( mouse );

		int fileMenu = glutCreateMenu( fileMenuFunc ) ;
		               glutAddMenuEntry( "ReLoad>Cube",    1 ) ;
		               glutAddMenuEntry( "SaveTo>NewCube", 2 ) ;
		       	       glutAddMenuEntry( "SaveTo>Slices",  3 ) ;
      		 
		int viewMenu = glutCreateMenu( viewMenuFunc ) ;
		               glutAddMenuEntry( "View>XY",  1 ) ;
		       	       glutAddMenuEntry( "View>XZ",  2 ) ;
		               glutAddMenuEntry( "View>YZ",  3 ) ;
		               glutAddMenuEntry( "View>XYZ", 4 ) ;
      		       
		int plotMenu = glutCreateMenu( plotMenuFunc ) ;
		               glutAddMenuEntry( "ByteHistogram",     1 ) ;
		               glutAddMenuEntry( "ShortHistogram",    2 ) ;
		               glutAddMenuEntry( "PlaneFocusCurve",   3 ) ;
		               glutAddMenuEntry( "PlaneMeanCurve",    4 ) ;
		               glutAddMenuEntry( "PlaneMaxCurve",     5 ) ;
		               glutAddMenuEntry( "PlotRefresh>On",    6 ) ;
		               glutAddMenuEntry( "PlotRefresh>Off",   7 ) ;
		               glutAddMenuEntry( "WriteTextFile>On",  8 ) ;
		               glutAddMenuEntry( "WriteTextFile>Off", 9 ) ;
      		        
		int operateMenu = glutCreateMenu( operateMenuFunc ) ;
		                  glutAddMenuEntry( "Contrast",       1 ) ;
		                  glutAddMenuEntry( "Padding",        2 ) ;
		                  glutAddMenuEntry( "Cropping",       3 ) ;
		                  glutAddMenuEntry( "Substract",      4 ) ;
		                  glutAddMenuEntry( "PSFshift",       5 ) ;
		                  glutAddMenuEntry( "lnTransform",    6 ) ;
		                  glutAddMenuEntry( "SpectraON/OFF",  7 ) ;
		                  glutAddMenuEntry( "LowPassFilter",  8 ) ;
      			 
		int controlMenu = glutCreateMenu( controlMenuFunc ) ;
		                  glutAddMenuEntry( "Xaxis>On/Off",   1 ) ;
		                  glutAddMenuEntry( "Yaxis>On/Off",   2 ) ;
		                  glutAddMenuEntry( "Zaxis>On/Off",   3 ) ;
		                  glutAddMenuEntry( "Xcross>On/Off",  4 ) ;
		                  glutAddMenuEntry( "Ycross>On/Off",  5 ) ;
		                  glutAddMenuEntry( "Zcross>On/Off",  6 ) ;
		                  glutAddMenuEntry( "XYlabel>On/Off", 7 ) ;
		                  glutAddMenuEntry( "XZlabel>On/Off", 8 ) ;
		                  glutAddMenuEntry( "YZlabel>On/Off", 9 ) ;
      
		glutCreateMenu( mainMenuFunc );
		glutAddSubMenu(   "File",    fileMenu    ) ;
		glutAddSubMenu(   "View",    viewMenu    ) ;
		glutAddSubMenu(   "Plot",    plotMenu    ) ;
		glutAddSubMenu(   "Operate", operateMenu ) ;
		glutAddSubMenu(   "Control", controlMenu ) ;
		glutAddMenuEntry( "Help",    1 ) ;
		glutAddMenuEntry( "Info",    2 ) ;
		glutAddMenuEntry( "Exit",    3 ) ;
		glutAttachMenu( GLUT_RIGHT_BUTTON ) ;
  		       
		glutMainLoop() ;
	}
	
	return 0 ;
}
