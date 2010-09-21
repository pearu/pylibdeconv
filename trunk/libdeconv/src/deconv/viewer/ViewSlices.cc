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
 * Filename:  ViewSlices.cc
 */


/*
 *	ViewSlices is designed for viewing and handling a series of 2-D sectioning images.
 *
 *	One input argument is required by the program: <images_setname>
 *	For example : $./ViewSlices /data/Seq0/Image
 */


#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <GL/glut.h>
#include "libdeconv/CCube.h"
#include "libdeconv/FFTW3fft.h"
#include "gnuplot_i.h"
#include "ViewerUtils.h"


#define WINDOW_UPPER_LEFT_X 250
#define WINDOW_UPPER_LEFT_Y 250


static int MainWindowID = 0 ;
static unsigned char* pixels ;

static gnuplot_ctrl * PlotWindow ;
static bool plot_refresh_on = true ;
static bool plot_txtfile_on = false ;

static double IntensityWindow ;
static double IntensityCenter ;
static int    MaxIntensity = 0 ;

static ShortSlice  image ;
static ShortSlice  frame ;

static char CurImageName[256] ;
static char CurFrameName[256] ;
static std::string images_setname ;
static std::string images_dirname ;
static std::vector< std::string > image_name ;

static int  CurImageID   = 0 ;
static int  ImageWidth   = 0 ;
static int  ImageHeight  = 0 ;
static int  ImageSize    = 0 ;
static int  FrameCenterX = 0 ;
static int  FrameCenterY = 0 ;
static int  FrameLeftX   = -1 ;
static int  FrameRightX  = -1 ;
static int  FrameTopY    = -1 ;
static int  FrameBottomY = -1 ;
static bool DisplayFrameWindow = false ;
  
static int  NumFastedImages = 1 ;


std::string usage = "$./ViewSlices <images_setname>\n" ;


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

	double* index = new double[ndata] ;
//	DoubleArray index( new double[ndata] ) ;
	for ( int i = 0 ; i < ndata ; i++ ) index[i] = (double) i ;

	if( plot_refresh_on ) gnuplot_resetplot( PlotWindow ) ;
	gnuplot_setstyle( PlotWindow, "lines" ) ;
	gnuplot_plot_xy( PlotWindow, index, data, ndata, title ) ;
	delete [] index;
}


void print_help( void )
{
	std::cout << "\n\n" ;
	std::cout << " ViewSlices Menu (see instruction in README).\n" ;
	std::cout << " ViewSlices Functional Keys: \n" ;
	std::cout << " ,           - fast   image(s).\n" ;
	std::cout << " .           - reward image(s).\n" ;
	std::cout << " f           - turn on/off frame window on viewer screen.\n" ;
	std::cout << " h           - print help.\n" ;
	std::cout << " ESC         - quit.\n" ;  
	std::cout << "\n\n" ;
}



void print_info( void )
{  	
	std::cout << "\n\n" ;
	std::cout << " Images loaded from directory  >>> " << images_dirname << "\n" ;
	std::cout << " Number of all searched images >>> " << image_name.size() << "\n" ;
	std::cout << " Image size : width x height   >>> " << ImageWidth << " x " << ImageHeight << "\n" ;
	std::cout << " Image type : max  intensity   >>> " << MaxIntensity << "\n" ;
	std::cout << " Number of fasted images once  >>> " << NumFastedImages << "\n" ;
	
	std::cout << " Intensity Contrasting Window  >>> " << IntensityWindow << "\n" ;
	std::cout << " Intensity Contrasting Center  >>> " << IntensityCenter << "\n" ;
	
	std::cout << " Image : " << CurImageName << ".\n" ;
	std::cout << "       max  intensity >>> " << image.slicemaxval() << "\n" ;
	std::cout << "       min  intensity >>> " << image.sliceminval() << "\n" ;
	std::cout << "       mean intensity >>> " << image.slicemeanval() << "\n" ;
	
	if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
	{
		std::cout << " SubImage : " << CurFrameName 
		          << " centered at Image(" << FrameCenterX << "," << FrameCenterY << ").\n" ;
		std::cout << "       width x height >>> " << frame.width() << " x " << frame.height() << "\n" ;
		std::cout << "       max  intensity >>> " << frame.slicemaxval() << "\n" ;
		std::cout << "       min  intensity >>> " << frame.sliceminval() << "\n" ;
		std::cout << "       mean intensity >>> " << frame.slicemeanval() << "\n" ;
	}
	std::cout << "\n\n" ;
}



int order_images( int& StartID, int& NumCubedImages, int& SkippedImages )
{
	std::cout << " Now order images to the cube, the ImageID ranges from 0 to " << image_name.size()-1 << ".\n" ; 
	std::cout << " Example of ordering images to the cube : \n" ; 
	std::cout << "        the first ImageID to   the cube = 0\n" ; 
	std::cout << "        number of images  to   the cube = 5\n" ; 
	std::cout << "        number of images  to be shipped = 2\n" ;
	std::cout << "        the cube will be made of images : image0, image3, image6, image9, image12\n" ;
	std::cout << " Input the first ImageID into the cube <<<< " ;
	std::cin  >> StartID ;
	std::cout << " Input number of images  into the cube <<<< " ;
	std::cin  >> NumCubedImages ;
	std::cout << " Input number of images  to be skipped <<<< " ;
	std::cin  >> SkippedImages ;
	
	if( StartID >= 0 && ( StartID+(NumCubedImages-1)*(SkippedImages+1) ) < (int) image_name.size() )
	{
		return 1 ;
	}
	else	return 0 ;
}



void write_images_to_cube( void )
{
	int        sliceID, StartID, NumCubedImages, SkippedImages ;
	char       cubename[256] ;
	ShortSlice sslice ;
	ByteSlice  bslice ;
	
	if( order_images( StartID, NumCubedImages, SkippedImages ) )
	{	
		std::cout << " Input the cube name without suffix <<<< " ;
		std::cin  >> cubename ;
		if( MaxIntensity == 65535 )
		{
			ShortCube  scube( ImageWidth, ImageHeight, NumCubedImages ) ;
			for( int i = 0 ; i < NumCubedImages ; i++ )
			{
				sliceID = StartID + i * ( SkippedImages + 1 ) ;
				sslice.read( image_name[sliceID].c_str() ) ;
				scube.setslice( i, XY_SLICE, sslice ) ;
				std::cout << " Image: " << image_name[sliceID] << " -> Cube: " << cubename << ".\n" ;
			}
			scube.write( cubename ) ;
		}
		else
		{
			ByteCube  bcube( ImageWidth, ImageHeight, NumCubedImages ) ;
			for( int i = 0 ; i < NumCubedImages ; i++ )
			{
				sliceID = StartID + i * ( SkippedImages + 1 ) ;
				bslice.read( image_name[sliceID].c_str() ) ;
				bcube.setslice( i, XY_SLICE, bslice ) ;
				std::cout << " Image: " << image_name[sliceID] << " -> Cube: " << cubename << ".\n" ;
			}
			bcube.write( cubename ) ;
		}
	}
	else	std::cout << " Failed to write images to the cube for an incorrect order.\n " ;
}



void write_frames_to_cube( void )
{
	int        sliceID, StartID, NumCubedImages, SkippedImages ;
	char       cubename[256] ;
	ShortSlice sslice ;
	ByteSlice  bslice ;

	if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
	{
		if( order_images( StartID, NumCubedImages, SkippedImages ) ) 
		{
			std::cout << " Input the cube name without suffix <<<< " ;
			std::cin  >> cubename ;
			if( MaxIntensity == 65535 ) 
			{
				ShortCube  scube( frame.width(), frame.height(), NumCubedImages ) ;
				for( int i = 0 ; i < NumCubedImages ; i++ )
				{
					sliceID = StartID + i * ( SkippedImages + 1 ) ;
					sslice.read( image_name[sliceID].c_str() ) ;
					sslice.crop( frame.width(), frame.height(), FrameCenterX, FrameCenterY ) ;
					scube.setslice( i, XY_SLICE, sslice ) ;
					std::cout << " SubImage: " << image_name[sliceID] << " -> Cube: " << cubename << ".\n" ;
				}
				scube.write( cubename ) ;
			}
			else
			{
				ByteCube  bcube( frame.width(), frame.height(), NumCubedImages ) ;
				for( int i = 0 ; i < NumCubedImages ; i++ )
				{
					sliceID = StartID + i * ( SkippedImages + 1 ) ;
					bslice.read( image_name[sliceID].c_str() ) ;
					bslice.crop( frame.width(), frame.height(), FrameCenterX, FrameCenterY ) ;
					bcube.setslice( i, XY_SLICE, bslice ) ;
					std::cout << " SubImage: " << image_name[sliceID] << " -> Cube: " << cubename << ".\n" ;
				}
				bcube.write( cubename ) ;
			}
		}
		else	std::cout << " Failed to write subimages to the cube for an incorrect order.\n" ;
	}		
	else	std::cout << " Failed to write sub_images to the cube for SubImage has not been set up yet.\n" ;
}



void destroy( void )
{
  	gnuplot_close( PlotWindow ) ;
  	glutDestroyWindow( MainWindowID ) ;
  	exit( 0 ) ;
}



void initImage( void )
{
	image.read( image_name[ CurImageID ].c_str() ) ;
	std::string::size_type pos = image_name[CurImageID].find_last_of( '/' ) ; 
	std::string s = image_name[CurImageID].substr( pos+1, image_name[CurImageID].find_last_of( '.' )-pos-1 ) ;
	sprintf( CurImageName, "%s", s.c_str() ) ;
	
	if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
	{
		frame = image ;
		frame.crop( FrameRightX-FrameLeftX+1, FrameBottomY-FrameTopY+1, FrameCenterX, FrameCenterY ) ;
		sprintf( CurFrameName, "%s_%d_%d_%d_%d", CurImageName, FrameLeftX, FrameRightX, FrameTopY, FrameBottomY ) ;
	}
}



void initDisplay( void )
{
	unsigned char pixel ;

	for( int i = 0 ; i < ImageSize ; i++ )
	{
		pixel = byte_contrast( (double) (image.data())[i], IntensityWindow, IntensityCenter ) ;
		pixels[ i*3     ] = pixel ;
		pixels[ i*3 + 1 ] = pixel ;
		pixels[ i*3 + 2 ] = pixel ;
	}
}


		
void display( void )
{
	glClear( GL_COLOR_BUFFER_BIT ) ;
	glRasterPos2i( 0, ImageHeight-1 ) ;
	glPixelZoom( 1.0, -1.0 ) ;
	glDrawPixels( ImageWidth, ImageHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels ) ;
	if( DisplayFrameWindow )
	{
		glColor3f( 1.0, 0.0, 0.0 ) ;
    		glBegin( GL_LINES ) ;
    		  glVertex2i( FrameLeftX,  ImageHeight-FrameTopY ) ;
    		  glVertex2i( FrameRightX, ImageHeight-FrameTopY ) ;
    		glEnd() ;
    		glBegin( GL_LINES ) ;
    		  glVertex2i( FrameLeftX,  ImageHeight-FrameBottomY ) ;
    		  glVertex2i( FrameRightX, ImageHeight-FrameBottomY ) ;
    		glEnd() ;
    		glBegin( GL_LINES ) ;
    		  glVertex2i( FrameLeftX, ImageHeight-FrameTopY ) ;
    		  glVertex2i( FrameLeftX, ImageHeight-FrameBottomY ) ;
    		glEnd() ;
    		glBegin( GL_LINES ) ;
    		  glVertex2i( FrameRightX, ImageHeight-FrameTopY ) ;
    		  glVertex2i( FrameRightX, ImageHeight-FrameBottomY ) ;
    		glEnd() ;
	}
	glFlush() ;
}



void reshape( int w, int h )
{
  	glViewport( 0, 0, (GLsizei)w, (GLsizei)h ) ;
  	glMatrixMode( GL_PROJECTION ) ;
  	glLoadIdentity() ;
	gluOrtho2D( 0.0, w, 0.0, h ) ;
	glMatrixMode( GL_MODELVIEW ) ;
  	glLoadIdentity() ;
	glutReshapeWindow( ImageWidth, ImageHeight ) ;
}



void keyboard(unsigned char key, int x, int y)
{
	switch( key ) 
	{
		case 44: // ',' - reward
			if( CurImageID < NumFastedImages ) 
			{
				CurImageID = CurImageID - NumFastedImages + image_name.size()  ;
			}
			else	
			{
				CurImageID = CurImageID - NumFastedImages ;
			}
			initImage() ;
			initDisplay() ;
			printf(" Image %4d : %s", CurImageID, CurImageName);
//			std::cout << " Image : " << CurImageName << ".\n" ;
			glutPostRedisplay() ;
			break ;
			
		case 46: // '.' - fast
			if( (CurImageID + NumFastedImages) > ((int) image_name.size() - 1) )
			{
				CurImageID = CurImageID + NumFastedImages - image_name.size() ;
			}
			else
			{
				CurImageID = CurImageID + NumFastedImages ;
			}
			initImage() ;
			initDisplay() ;
			printf(" Image %4d : %s", CurImageID, CurImageName);
//			std::cout << " Image : " << CurImageName << ".\n" ;
			glutPostRedisplay() ;
			break ;
		
		case 102: // 'f' - frame widnow display
			if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
			{
				DisplayFrameWindow = !DisplayFrameWindow ;
				glutPostRedisplay() ;
			}
			else	std::cout << " SubImage has not been set up yet.\n" ;
			break ;
    			
		case 104: // h - help
			print_help() ;
			break ;

		case 27: // ESC - Quit
			destroy() ;
			break ;

		default:
			break ;
	}
}



void mouse(int button, int state, int x, int y)
{
	unsigned char curpixel[3] ;

	if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) 
	{
		glReadPixels( x, y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, curpixel ) ;
		printf( "Image(%4d, %4d) = %d\n", x, y, image(x, y) ) ;
	}
}



void fileMenuFunc( int choice )
{
	switch( choice ) 
	{
		case 1: 
			write_images_to_cube() ;
			break ;

		case 2: 
			write_frames_to_cube() ;
			break ;

		default:
			break ;
	}
}



void plotMenuFunc( int choice )
{
	int    maxindex = 0 ;
	double  maxdata = 0.0 ;
	double  sumdata = 0.0 ;
	int       ndata = MaxIntensity+1 ;
	int remove_peak = 0 ; 
	char title[256] ;
    	    	
	switch( choice ) 
	{
		case 1:
		{
			double* data = new double[ ndata ] ;
			for( int i = 0 ; i < ndata ; i++ )        data[i] = 0.0 ;    			
			for( int i = 0 ; i < image.size() ; i++ ) data[ (int) (image.data())[i] ]++ ;
			for( int i = 0 ; i < ndata ; i++ )
			{
				sumdata += data[i] ;
				if( data[i] >= maxdata ) 
				{ 
					maxdata = data[i] ; 
					maxindex = i ; 
				}
			}
			std::cout << " Histogram peak found at " << maxindex << " : " << maxdata << "/" << sumdata << "\n" ;    			
			std::cout << " Indicate whether removing(1) the peak or not(0) to plot data <<<< " ;
			std::cin  >> remove_peak ;
			if( remove_peak ) data[ maxindex ] = 0.0 ;
			gnu_plot_data( data, ndata, CurImageName ) ;
			delete [] data ;
			break ;
		}

		case 2:
		{
			if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
			{
				double* data = new double[ ndata ] ;
				for( int i = 0 ; i < ndata ; i++ )        data[i] = 0.0 ; 
				for( int i = 0 ; i < frame.size() ; i++ ) data[ (int) (frame.data())[i] ]++ ;
				for( int i = 0 ; i < ndata ; i++ )
				{
					sumdata += data[i] ;
					if( data[i] >= maxdata ) 
					{ 
						maxdata = data[i] ; 
						maxindex = i ; 
					}
				}
				std::cout << " Histogram peak found at " << maxindex << " : " << maxdata << "/" << sumdata << "\n" ;    			
				std::cout << " Indicate whether removing(1) the peak or not(0) to plot data <<<< " ;
				std::cin  >> remove_peak ;
				if( remove_peak ) data[ maxindex ] = 0.0 ;
				gnu_plot_data( data, ndata, CurImageName ) ;
				delete [] data ;
			}
			else 	std::cout << " SubImage has not been set up yet.\n" ;
			
			break ;
		}
		
		case 3:
		{
			DoubleSlice dslice( ImageWidth, ImageHeight ) ;
			double* data = new double[ image_name.size() ] ;
			for( int k = 0 ; k < ((int) image_name.size()) ; k++ ) 
			{
				data[k] = 0.0 ;
				dslice.read( image_name[k].c_str() ) ;
				for( int j = 0 ; j < dslice.height() ; j++ )
					for( int i = 0 ; i < dslice.width() ; i++ )
					{
						if( i < dslice.width()-1 ) 
						{
							data[k] += ( dslice(i,j) * ( dslice(i,j)-dslice(i+1,j) ) ) ;
						}
					}
				if( data[k] >= maxdata ) 
				{ 
					maxdata = data[k] ; 
					maxindex = k ; 
				}
			}
			std::cout << "Sectioning is Focused at "<< image_name[ maxindex ] << " : " << maxdata << "\n" ;
			sprintf( title, "Focus along Sectioning Images" ) ;
			gnu_plot_data( data, image_name.size(), title ) ;
			delete [] data ;
			break;
		}
		
		case 4:
		{
			if( FrameBottomY > FrameTopY && FrameRightX > FrameLeftX )
			{
				DoubleSlice dslice( ImageWidth, ImageHeight ) ;
				double* data = new double[ image_name.size() ] ;
				for( int k = 0 ; k < ((int) image_name.size()) ; k++ ) 
				{
					data[k] = 0.0 ;
					dslice.read( image_name[k].c_str() ) ;
					dslice.crop( FrameRightX-FrameLeftX+1, FrameBottomY-FrameTopY+1, FrameCenterX, FrameCenterY ) ;
					for( int j = 0 ; j < dslice.height() ; j++ )
						for( int i = 0 ; i < dslice.width() ; i++ )
						{
							if( i < dslice.width()-1 ) 
							{
								data[k] += ( dslice(i,j) * ( dslice(i,j)-dslice(i+1,j) ) ) ;
							}
						}
					if( data[k] >= maxdata ) 
					{ 
						maxdata = data[k] ; 
						maxindex = k ; 
					}
				}
				std::cout << "SubImage Sectioning is Focused at "<< image_name[ maxindex ] << " : " << maxdata << "\n" ;
				sprintf( title, "Focus_along_Sectioning_SubImages_%d_%d_%d_%d", FrameLeftX, FrameRightX, FrameTopY, FrameBottomY ) ;
				gnu_plot_data( data, image_name.size(), title ) ;
				delete [] data ;
			}
			else 	std::cout << " SubImage has not been set up yet.\n" ;
			break ;
		}	
    			
		case 5:
			plot_refresh_on = true ;
			break ;
    			
		case 6:
			plot_refresh_on = false ;
			break ;
    			
		case 7:
			plot_txtfile_on = true ;
			break ;
    		
		case 8:
			plot_txtfile_on = false ;
			break ;
    	
		default:
			break;
	}
}


void operateMenuFunc( int choice )
{
	double value ;
	int    x, y, w, h ;

	switch( choice ) 
	{
		case 1:
			std::cout << " Input the new Intensity_Window_Width  <<<< " ;
			std::cin  >> value ;
			if( value <= 0.0 )
			{
				std::cout << " Wrong Input : Intensity_Window_Width must be positive.\n" ;
			}
			else
			{
				IntensityWindow = value ;
				std::cout << " Input the new Intensity_Window_Center <<<< " ;
				std::cin  >> IntensityCenter ;
				initDisplay() ;
				glutPostRedisplay() ;
			}
			break ;

    		case 2:
			std::cout << " The ImageID ranges from 0 to " << image_name.size()-1 << ".\n" ;
			std::cout << " Input the ImageID to display <<<< " ;
			std::cin  >> x ;
			if( x >= 0 && x < ((int) image_name.size()) )
			{ 
				CurImageID = x ;
				initImage() ;
				initDisplay() ;
				glutPostRedisplay() ;
			}
			else	std::cout << " Wrong Input : inputed ImageID is out of range.\n" ;		
			break ;     			

		case 3:
			std::cout << " Input number of images to be fasted/rewarded once <<<< " ;
			std::cin  >> x ;
			if( x > 0 )
			{ 
				NumFastedImages = x ;
			}
			else	std::cout << " Wrong Input : number must be positive.\n" ;
			break ;

		case 4:
		{
			std::cout << " Image size (width x height) = " << ImageWidth << "x" << ImageHeight << "\n" ;
			std::cout << " Input SubImage centered at image X-coordinate <<<< " ;
			std::cin  >> x ;
			std::cout << " Input SubImage centered at image Y-coordinate <<<< " ;
			std::cin  >> y ;
			std::cout << " Input SubImage width  (must be even) <<<< " ;
			std::cin  >> w ;
			std::cout << " Input SubImage height (must be even) <<<< " ;
			std::cin  >> h ;
			if( w%2 == 0 && h%2 == 0 )
			{
				ShortSlice sslice = image ;
				if( !sslice.crop( w, h, x, y ) )
				{
					FrameCenterX = x ;
					FrameCenterY = y ;
					FrameLeftX   = FrameCenterX - w/2 ; 
					FrameRightX  = FrameCenterX + w/2 - 1 ;
					FrameTopY    = FrameCenterY - h/2 ; 
					FrameBottomY = FrameCenterY + h/2 - 1 ;
					sprintf( CurFrameName, "%s_%d_%d_%d_%d", CurImageName, 
					         FrameLeftX, FrameRightX, FrameTopY, FrameBottomY ) ;
					frame = sslice ;
					std::cout << " SubImage: " << CurFrameName << " is set and centered at Image(" 
					          << FrameCenterX << "," << FrameCenterY << ").\n" ;
				}
				else	std::cout << " Wrong input : SubImage is out of image size.\n" ;
			}
			else	std::cout << " Wrong input : SubImage width or height must be even.\n" ;
			break ;
		}

		default:
			break;
	}
}



void mainMenuFunc( int choice) 
{
	switch( choice ) 
	{

		case 1:
			print_info() ;
			break ;

		case 2:
			print_help() ;
			break ;

		case 3:
			destroy() ;
			break ;

		default:
			break ;
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
		images_setname = std::string( argv[1] ) ;
		images_dirname = images_setname.substr( 0, images_setname.find_last_of( '/' ) ) ;
	}
 
	struct dirent * dp = NULL ;
	DIR * dirp = opendir( images_dirname.c_str() ) ;
	while( dirp && (dp = readdir(dirp)) != NULL )
	{
		std::string filename = images_dirname + '/' + std::string( dp->d_name ) ;
		if( filename.size() >= images_setname.size() )
		{
			if( filename.compare( 0, images_setname.size(), images_setname ) == 0 ) 
			{
				image_name.push_back( filename ) ;
			}
		}
	}
	
	if( image_name.size() > 0 )
	{
		std::cout << " There are " << image_name.size() << " images found. Now checking them ... " ;
		MaxIntensity = image.read( image_name[0].c_str() ) ;
		ImageWidth = image.width() ;
		ImageHeight = image.height() ;
		ImageSize = ImageWidth * ImageHeight ;
		for( int i = 1 ; i < (int) image_name.size() ; i++ )
		{
			if( image.read( image_name[i].c_str() ) != MaxIntensity || 
			    image.width() != ImageWidth || image.height() != ImageHeight )
			{
				std::cout << " Program is terminated for loading an unexpected image : " 
				          << image_name[i] << "\n" ;
				std::cout << " All of images must have same size and same type (byte or short).\n" ; 
				return 1 ;
			}
		}
		std::cout << " done.\n" ;		

		CurImageID = (int) ((floor) (image_name.size() / 2)) ;
		IntensityWindow = (double) MaxIntensity ;
		IntensityCenter = IntensityWindow / 2.0 ;
		pixels = new unsigned char [ ImageSize * 3 ] ;
		PlotWindow = gnuplot_init() ;
		initImage() ;
		initDisplay() ;

		glutInit( &argc, argv ) ;
		glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB ) ;
		glutInitWindowSize( ImageWidth, ImageHeight ) ;
		glutInitWindowPosition( WINDOW_UPPER_LEFT_X, WINDOW_UPPER_LEFT_Y ) ;
		glClearColor( 0.0, 0.0, 0.0, 0.0 ) ;
		glClear( GL_COLOR_BUFFER_BIT ) ;
		glShadeModel( GL_FLAT ) ;
		glPixelStorei( GL_UNPACK_ALIGNMENT, 1 ) ;
		MainWindowID = glutCreateWindow( images_setname.c_str() ) ;
		glutDisplayFunc( display ) ;
		glutReshapeFunc( reshape ) ;
		glutKeyboardFunc( keyboard ) ;
		glutMouseFunc( mouse ) ;

		int fileMenu = glutCreateMenu( fileMenuFunc ) ;
		               glutAddMenuEntry( "Images>Cube",    1 ) ;
		               glutAddMenuEntry( "SubImages>Cube", 2 ) ;
 
		int plotMenu = glutCreateMenu( plotMenuFunc ) ;
		               glutAddMenuEntry( "Image>Histogram",    1 ) ;
		               glutAddMenuEntry( "SubImage>Histogram", 2 ) ;
		               glutAddMenuEntry( "Image>Focus",        3 ) ;
		               glutAddMenuEntry( "SubImage>Focus",     4 ) ;
		               glutAddMenuEntry( "PlotRefresh>On",     5 ) ;
		               glutAddMenuEntry( "PlotRefresh>Off",    6 ) ;
		               glutAddMenuEntry( "WriteTextFile>On",   7 ) ;
		               glutAddMenuEntry( "WriteTextFile>Off",  8 ) ;

  		int operateMenu = glutCreateMenu( operateMenuFunc ) ;
		                  glutAddMenuEntry( "Set>Contrast", 1 ) ;
		                  glutAddMenuEntry( "Set>ImageID",  2 ) ;                  
		                  glutAddMenuEntry( "Set>Fast",     3 ) ; 
		                  glutAddMenuEntry( "Set>SubImage", 4 ) ;
                    
		glutCreateMenu( mainMenuFunc ) ;
		glutAddSubMenu( "File",    fileMenu ) ;
		glutAddSubMenu( "Plot",    plotMenu ) ;
		glutAddSubMenu( "Operate", operateMenu ) ;
		glutAddMenuEntry( "Info", 1 ) ;
		glutAddMenuEntry( "Help", 2 ) ;
		glutAddMenuEntry( "Exit", 3 ) ;    
		glutAttachMenu( GLUT_RIGHT_BUTTON ) ;

		glutMainLoop() ;
	}
	else
	{
		std::cout << " No " << images_setname << "* images found.\n" ;  
		return 1 ;
	}

	delete [] pixels ;
	return 0 ;
}
