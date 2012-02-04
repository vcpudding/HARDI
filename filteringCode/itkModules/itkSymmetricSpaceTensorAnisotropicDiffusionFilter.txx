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
#ifndef __itkSymmetricSpaceTensorAnisotropicDiffusionFilter_txx_
#define __itkSymmetricSpaceTensorAnisotropicDiffusionFilter_txx_

#include "itkSymmetricSpaceTensorAnisotropicDiffusionFilter.h"

//helper method
template<class T>
extern bool IsValidTensor(T &tensor);



namespace itk{

template <class TInputImage, class TOutputImage>
void
SymmetricSpaceTensorAnisotropicDiffusionFilter<TInputImage, TOutputImage>
::ThreadedApplyUpdate(TimeStepType dt, const ThreadRegionType &regionToProcess,
                      int threadId)
{
  std::cout<<"Inside Threaded Apply Update"<<std::endl;	
  ImageRegionIterator<UpdateBufferType> u(this->GetUpdateBuffer(),regionToProcess);
  ImageRegionIterator<OutputImageType>  o(this->GetOutput(), regionToProcess);

  u = u.Begin();
  o = o.Begin();

  TensorType refTensor,updatedTensor,deltaTensor;
  PixelType updatedValue,originalValue,deltaValue;
  SymmetricSpaceTensorGeometry<typename TOutputImage::PixelType::ValueType, 3>  ssTensorGeometry;

  long updatedCount=0;
  while ( !u.IsAtEnd() )
  {
	  originalValue=o.Value();
	  deltaValue= static_cast<PixelType>(u.Value() * dt);
	  for(int i=0;i<6;i++)
	  {
		  refTensor[i]=originalValue[i];
		  deltaTensor[i]=deltaValue[i];
	  }
	  if(IsValidTensor(refTensor)) {
		  updatedTensor=ssTensorGeometry.ExpMap(refTensor,deltaTensor);
		  for(int i=0;i<6;i++)
			  updatedValue[i]=updatedTensor[i];
	  }

	  if(IsValidTensor(refTensor)){
		  o.Value()=updatedValue ;  // no adaptor support here
		  updatedCount++;
	  }
	  ++o;++u;
  }
  std::cout<<"Updated count="<<updatedCount<<std::endl;
}

}//end namespace itk

#endif
