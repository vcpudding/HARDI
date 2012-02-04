/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSymmetricSpaceTensorAnisotropicDiffusionFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/12/09 12:15:42 $
  Version:   $Revision: 1.28 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianFilter_txx_
#define __itkGaussianFilter_txx_
#include<iomanip>
#include<fstream>
#include<ios>
namespace itk{

template <class TInputImage, class TOutputImage>
void
GaussianFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      int threadId)
{

  std::cout<<"Inside Gaussian Threaded Apply Update"<<std::endl;	
  //delta is the smoothing prior
  //u0 is the noisy image being filtered
  //u is the filtered image(as of current iteration)
  //We add an attachment term to the delta 
  
  ImageRegionIterator<UpdateBufferType> delta(this->GetUpdateBuffer(),regionToProcess);
  ImageRegionIterator<OutputImageType>  u0(noisyImage,regionToProcess);
  ImageRegionIterator<OutputImageType>  u(this->GetOutput(), regionToProcess);
  
  delta = delta.Begin();
  u0=u0.Begin();
  u = u.Begin();

  double _u0,_u,_delta,_bias,_netDelta=0;

   while ( !delta.IsAtEnd() )
   {
/* 	    std::cout<<"u Before="<<u.Value()<<std::endl; */
/* 	    std::cout<<"u0 Before="<<u0.Value()<<std::endl; */
/* 	    std::cout<<"Delta ="<<delta.Value()<<std::endl; */
	   for(int i=0;i<(u.Value()).Size();i++)
	   {
		   _u=static_cast<double>((u.Value())[i]);
		   _u0=static_cast<double>((u0.Value())[i]); 
		   _delta=static_cast<double>((delta.Value())[i] );
		   
		   //use a gaussain correction term
		   _bias=at->GetAttachmentTerm(_u0,_u);
		   _netDelta= ((_delta + _bias)*dt);
		   if(u.Value()[i])
		   (u.Value())[i]+=_netDelta ;  // no adaptor support here
		  // std::cout<<"_u="<<_u<<",_u0="<<_u0<<std::endl;
		  // std::cout<<"delta="<<_delta<<",bias="<<_bias<<std::endl;
	   }
/* 	     std::cout<<"u After="<<u.Value()<<std::endl; */
/* 	     std::cout<<"u0 After="<<u0.Value()<<std::endl; */
	   ++u0;
	   ++u;
	   ++delta;
   }

  std::cout<<"After Gaussian Apply Update"<<std::endl;

}
}//end namespace itk

#endif
