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
 * Filename:  CGdeconvolver.h
 */
 

#ifndef CGDECONVOLVER_H
#define CGDECONVOLVER_H


#include "LWCGdeconvolver.h"


/*
 *	======================================================================================================
 *	CGdeconvolver is developped based on the Maximum Likelihood-Conjugate Gradient Iterative Deconvolution 
 *	(ML-CGID) method, and is used for the 3-D deconvolution applied to fluorescence microscopy imaging. 
 *	======================================================================================================	
 *
 *
 *		----------------------------------------------------------
 *		CG Intensity Regularization : <_IsApplyIR>, <_CGIRpenalty>
 *		----------------------------------------------------------
 *		Intensity regularization will be applied iteratively if <_IsApplyIR> is set to be true. 
 *		The intensity penalty <_CGIRpenalty> will be calculated in CGdeconvolver::run(). 
 *	
 *
 *		----------------------------------------------------------------------
 *		run() : <image>, <psf>, <object>, <SpacialSupport>, <FrequencySupport>
 *		----------------------------------------------------------------------
 *
 *		*****************************
 *		The Input 3-D Image : <image>  
 *		*****************************
 *		The input cubic image to be deconvolved, with the dimensions of DimX, DimY and DimZ 
 *		where DimX is its fastest varying dimension and DimZ is its slowest varying dimension, 
 *		must have a black (0) background. Each dimension of the cubic image must be a power of 2. 
 *		The image data must be stored in an one-dimensional array : <image> as "x+y*DimX+z*DimX*DimY".
 *		
 *		Warnning: the <image> array will be rewritten in run(). 
 *
 *		*************************
 *		The Input 3-D PSF : <psf>
 *		*************************
 *		The 3-D PSF must have the same dimensions as the input image and its data must be stored 
 *		in an one-dimensional array : <psf> in a special manner for the deconvolution.
 *		
 *		Generally, a 3-D PSF has an hour-glass shape, however the hour-glass shaped PSF has to be shifted 
 *		in a special manner for deconvolution. This is termed as its deconvolution shape here. 
 *
 *		A 3-D PSF, generated from the Fluo3DPSF::create() or FluoRZPSF::get3Dpsf(), 
 *		has been shifted for deconvolution and can be directly taken into run().
 *		Make sure that the input PSF is stored correctly if you create it in your own way.
 *		Use functions "Operate->PSFshift and Operate->lnTransform" in the "ViewCube" viewer 
 *		to see how a 3-D PSF is shifted between the hour-glass and deconvolution shapes.  
 *
 *		Warnning: the <psf> array will be rewritten in run(). 
 *	
 *		************************************************************************
 *		The First Estimated Object and The Finally Deconvolved Object : <object>
 *		************************************************************************
 *		The first estimated object could be same as the input cubic image or a previously deconvolved object
 *		on the input cubic image. The latter case is used for you want to continue a deconvolution process on 
 *		a same cubic image from a previously deconvolved object on this image. The previously deconvolved 
 *		object does not have to be obtained in a same deconvolution routine as used this time.
 *
 *		<object> must be an one-dimensional array storing the the first estimated object data and also used 
 *		to store the finally deconvolved object data. Both are stored in a same manner as the input image.
 *
 *		**************************************
 *		The Spacail Support : <SpacialSupport>
 *		**************************************
 *		Spacial support is not applied by default, but can be applied on estimated objects iteratively
 *		by passing an one-dimensional unsigned char array <SpacialSupport>, whose dimensions are same
 *		as the input image, to run(). Each value in this array must be 0 or 1.
 *
 *		******************************************
 *		The Frequency Support : <FrequencySupport>
 *		******************************************
 *		Frequency support is not applied by default, but can be applied on the FT of the PSF by passing
 *		an one-dimensional unsigned char array <FrequencySupport>, whose dimensions are same as the PSF, 
 *		to run(). Each value in this array must be 0 or 1.
 */



/*
 *	CGdeconvolver working space in double floating precision
 */
typedef struct
{
	int size ;
	double * psf_re ;
	double * psf_im ;
	double * image_re ;
	double * image_im ;
	double * cg_re ;
	double * cg_im ;
	double * otf ;
	double * object0 ;
	unsigned char * sign ;
} CGdws ;


/*
 *	CGdeconvolver working space in single floating precision
 */
typedef struct
{
	int size ;
	float * psf_re ;
	float * psf_im ;
	float * image_re ;
	float * image_im ;
	float * cg_re ;
	float * cg_im ;
	float * otf ;
	float * object0 ;
	unsigned char * sign ;
} CGsws ;

               
class CGdeconvolver : public LWCGdeconvolver
{
	public:
//	typedef boost::shared_ptr< CGdeconvolver > Ptr ;	
		~CGdeconvolver() ;
		CGdeconvolver () ;
	
	/*
	 *	Set up the control flags and default parameters used for CGdeconvolver
	 *	Input:
	 *		IsApplyIR,     it is to indicate whether to apply intensity regularization or not;
	 *		               see the description for the CG intensity regularization above in this file.
	 *		IsApplyNorm,   it is the apply_normalization   indicator. (described in "deconvolver.h")
	 *		IsTrackLike,   it is the track_likelihood      indicator. (described in "deconvolver.h")
	 *		IsTrackMax,    it is the track_max_intensity   indicator. (described in "deconvolver.h")
	 *		IsCheckStatus, it is the check_program_running indicator. (described in "deconvolver.h")
	 *	Warning:
	 *		<_MaxRunIteration>       will be set to be its default value. (described in "deconvolver.h")
	 *		<_Criterion>             will be set to be its default value. (described in "deconvolver.h")
	 *		<_ConditioningIteration> will be set to be its default value. (described in "LWCGdeconvolver.h")
	 *		<_ConditioningValue>     will be set to be its default value. (described in "LWCGdeconvolver.h")
	 *		<_ConditioningTolerance> will be set to be its default value. (described in "LWCGdeconvolver.h")
	 */ 
	void    init( bool IsApplyIR = true,   bool IsApplyNorm = false, bool IsTrackLike = false, 
	              bool IsTrackMax = false, bool IsCheckStatus = true ) ;
        
        
	/*
	 *	Run CGdeconvolution in the double/single floating precision
	 *	Input:
	 *		DimX,             it is the fastest varying dimension of the image/psf; it must be power of 2.
	 *		DimY,             it is the middle          dimension of the image/psf; it must be power of 2.
	 *		DimZ,             it is the slowest varying dimension of the image/psf; it must be power of 2.
	 *		image,            it points to an one-dimensional DimX*DimY*DimZ array storing the image data.
	 *		psf,              it points to an one-dimensional DimX*DimY*DimZ array storing the psf data.
	 *		object,           it points to an one-dimensional DimX*DimY*DimZ array storing both 
	 *		                  the first estimated object data and the finally deconvolved object data. 
	 *		ws,               it points to the CGdeconvolver double/float working space. 
	 *		SpacialSuppport,  it points to an one-dimensional unsigned char DimX*DimY*DimZ
	 *		                  array storing the spacial support data and its default is NULL.
	 *		FrequencySupport, it points to an one-dimensional unsigned char DimX*DimY*DimZ
	 *		                  array storing the frequency support data and its default is NULL.
	 *	Throw:
	 *		throw an error if a given dimension is wrong.
	 */                   
	void    run( int DimX, int DimY, int DimZ, double * image, double * psf, double * object, CGdws & ws, 
	             unsigned char * SpacialSupport = NULL, unsigned char * FrequencySupport = NULL ) ;   
	void    run( int DimX, int DimY, int DimZ, float  * image, float  * psf, float  * object, CGsws & ws,
	             unsigned char * SpacialSupport = NULL, unsigned char * FrequencySupport = NULL ) ;
	
	/*
	 *	Export the profile of a CGdeconvolver to a text file
	 *	Input:
	 *		filename, it is the name of the text file to be written including suffix. 
	 *	Throw:
	 *		throw an error if fail.
	 */	  
	void    exportCG( const char * filename ) ;
	

	private:
	bool           	_ApplyIR ;
	double         	_CGIRpenalty ;	
	FFTW3_FFT* 		_FFTplanff ;
	FFTW3_FFT* 		_FFTplanbb ;
        
	void    _CGprintStatus( int stage ) ;
        
	void    _CGstartRun( int DimX, int DimY, int DimZ, CGdws & ws ) ;
	void    _CGstartRun( int DimX, int DimY, int DimZ, CGsws & ws ) ;
        
	void    _CGfinishRun( CGdws & ws ) ;
	void    _CGfinishRun( CGsws & ws ) ;
        
	void    _CGupdate1( double & gamma, double & alpha, double & beta, double * cgr, double * cgp, CGdws & ws ) ; 
	void    _CGupdate1( float  & gamma, float  & alpha, float  & beta, float  * cgr, float  * cgp, CGsws & ws ) ; 
        
	void    _CGupdate2( double & gamma, double & alpha, double & beta, double * cgr, double * cgp, CGdws & ws ) ; 
	void    _CGupdate2( float  & gamma, float  & alpha, float  & beta, float  * cgr, float  * cgp, CGsws & ws ) ;
  
	double  _CGrunRegularization( int size, double * img, double * otf ) ;
	double  _CGrunRegularization( int size, float  * img, float  * otf ) ;
} ;


#endif   /*   #include "CGdeconvolver.h"   */
