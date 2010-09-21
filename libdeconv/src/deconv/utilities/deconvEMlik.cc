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
 * Filename:  deconvEM.cc
 */
 

#include "libdeconv/EMdeconvolver.h"
#include "libdeconv/CCube.h"


/*
 *	deconvEM is designed for 3-D deconvolution using the EMdeconvolver and CCube classes.
 *
 *	Four required input arguments:
 *		<image_name> <psf_name> <first_estimated_object_name> <deconvolved_object_name_without suffix> <Acceleration>	
 *	Four dummy input arguments:
 *		[<max_allowed_deconvolved_iterations>] [<criterion_to_stop_deconvolution>] [IR_Iterations>] [IR_penalty>]
 *		
 *	The deconvolved object will be saved using CCube::save():
 *		The header file : <deconvolved_object_name_without suffix>.hdr
 *		The data   file : <deconvolved_object_name_without suffix>.f32 if <DataType> is float. 
 *				  <deconvolved_object_name_without suffix>.f64 if <DataType> is double.
 *
 *	The deconvolution profile will be saved in <deconvolved_object_name_without suffix>_plan.txt
 *
 *	The update data will be saved in <deconvolved_object_name_without suffix>_Update.txt
 *	
 *	If <Track_Likelihood> is defined to be true :
 *		The likelihood data will be saved in <deconvolved_object_name_without suffix>_Likelihood.txt
 *
 *	If <Track_Max_In_Object> is defined to be true:
 *		The max intensity data will be saved in <deconvolved_object_name_without suffix>_Likelihood.txt	
 */
 

#define DataType                    float   //deconvolution precision
#define WorkSpaceType               EMsws   //put EMdws if <DataType> is double.

#define Apply_Normalization         false   //flag to control whether to normalize the input image and the iterative estimated object
#define Track_Likelihood            true   //flag to control whether to track the likelihood value iteratively
#define Track_Max_In_Object         false   //flag to control whether to track max intensity of the iterarive estimated object
#define Check_Program_Run           true    //flag to control whether to check program running status in terminal


std::string usage = std::string("$deconvEM <image> <psf> <first_estimated_object> <deconvolved_object> <Acceleration>\n") + 
                    std::string("          [<max_allowed_deconvolved_iterations>] [<criterion_to_stop_deconvolution>]\n") +
                    std::string("          [<IR_Iterations>] [<IR_penalty>]\n") ;


int main( int argc, char ** argv )
{
	CCube<DataType>      image ;
	CCube<DataType>      psf ;
	CCube<DataType>      object ;
	EMdeconvolver*       em = new EMdeconvolver () ;
	WorkSpaceType        ws ;
	FILE*                fp ;
	char                 filename[256] ;
	double*              vtrack ;
	int                  vtrack_size ;
	bool                 acceleration ;
	

	if( argc >=6 && argc <= 10 )
	{
		image.read( argv[1] ) ;
		psf.read( argv[2] ) ;
		object.read( argv[3] ) ;
		
		if( atoi(argv[5]) == 0 ) acceleration = false ;
		else                     acceleration = true ;

		em->init( acceleration, Apply_Normalization, Track_Likelihood, Track_Max_In_Object, Check_Program_Run ) ; 
		
		if( argc > 6 ) em->setMaxRunIteration( (unsigned int)( atoi(argv[6]) ) ) ; 
		if( argc > 7 ) em->setCriterion( (double)( atof(argv[7]) ) ) ;   
		if( argc > 8 ) em->setEMIRiteration( (unsigned int)( atoi(argv[8]) ) ) ;
		if( argc > 9 ) em->setEMIRpenalty( (double)( atof(argv[9]) ) ) ;
		
		em->run( image.length(), image.width(), image.height(), image.data(), psf.data(), object.data(), ws ) ;
	
		sprintf( filename, "%s_plan.txt", argv[4] ) ;
		em->exportEM( filename ) ;
		
		vtrack = new double [ em->MaxRunIteration() ] ;

		vtrack_size = em->exportUpdateTrack( vtrack ) ;
		if( vtrack_size > 0 )
		{
			sprintf( filename, "%s_Update.txt", argv[4] ) ;
			fp = fopen( filename, "wt" ) ;
			for( int i = 0; i < vtrack_size; i++ ) fprintf( fp, "%4d %12.6e\n", i, vtrack[i] ) ;
			fclose( fp ) ;

			if( em->TrackMaxInObject() )
			{
				vtrack_size = em->exportObjectMaxTrack( vtrack ) ;
				sprintf( filename, "%s_MaxIntensity.txt", argv[4] ) ;
				fp = fopen( filename, "wt" ) ;
				for( int i = 0; i < vtrack_size; i++ ) fprintf( fp, "%4d %12.6e\n", i, vtrack[i] ) ;
				fclose( fp ) ;
			}
			else	std::cout << " MaxIntensity of deconvolved object wasn't tracked iteratively.\n" ;

			if( em->TrackLikelihood() )
			{
				vtrack_size = em->exportLikelihoodTrack( vtrack ) ;
				sprintf( filename, "%s_Likelihood.txt", argv[4] ) ;
				fp = fopen( filename, "wt" ) ;
				for( int i = 0; i < vtrack_size; i++ ) fprintf( fp, "%4d %12.6e\n", i, vtrack[i] ) ;
				fclose( fp ) ;
			}
			else	std::cout << " Likelihood wasn't tracked iteratively.\n" ;

			object.write( argv[4] ) ;
		}
		else	std::cout << " EMdeconvolver::run() didNOT go into the loop.\n" ;
		
		delete [] vtrack ;
		delete em ;
	}
	else
	{
		std::cout << usage ;
		return 1 ;
	}


	return 0 ;
}
		
	
		
