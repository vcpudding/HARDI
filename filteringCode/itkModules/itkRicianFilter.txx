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
#ifndef __itkRicianFilter_txx_
#define __itkRicianFilter_txx_
#include<iomanip>
#include<fstream>
#include<ios>
namespace itk{

template <class TInputImage, class TOutputImage>
void
RicianFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      int threadId)
{

  std::cout<<"Inside Rician Threaded Apply Update"<<std::endl;	
  //delta is the smoothing prior
  //u0 is the noisy image being filtered
  //u is the filtered image(as of current iteration)
  //We add an attachment term to the delta to account
  //for the rician bias
  ImageRegionIterator<UpdateBufferType> delta(this->GetUpdateBuffer(),regionToProcess);
  ImageRegionIterator<OutputImageType>  u0(noisyImage,regionToProcess);
  ImageRegionIterator<OutputImageType>  u(this->GetOutput(), regionToProcess);

  delta = delta.Begin();
  u0=u0.Begin();
  u = u.Begin();

  double _u0,_u,_delta,_bias,_netDelta;
  long count=0; 
  double delta_accum=0,bias_accum=0;
  long delta_counter=0;

  while ( !delta.IsAtEnd() )
  {
	  for(int i=0;i<(u.Value()).Size();i++)
	  {
		  _u=static_cast<double>((u.Value())[i]);
		  _u0=static_cast<double>((u0.Value())[i]); 
		  _delta=static_cast<double>((delta.Value())[i] );
		  delta_accum+=(_delta*_delta);
		  delta_counter++;
		  double sigma=at->Getsigma();

		  _bias=at->GetAttachmentTerm(_u0,_u); 
		  bias_accum+=_bias; 
		  _netDelta= (_delta + _bias)*dt;

		  if((u.Value())[i])
			  (u.Value())[i]+=_netDelta ;  // no adaptor support here
		  //std::cout<<"_u="<<_u<<",_u0="<<_u0<<std::endl;
		  //<"_u0<<",delta="<<_delta<<",bias="<<_bias<<std::endl;
	  }
	  //	  std::cout<<"u After="<<u.Value()<<std::endl;
	  //	 std::cout<<"u0 After="<<u0.Value()<<std::endl;

	  ++u0;
	  ++u;
	  ++delta;
  }
}
}//end namespace itk

#endif
