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


  long invalid=0;
    long count=0;
while ( !u.IsAtEnd() )
    {
    originalValue=o.Value();
/*     if(isnan(dt)) */
/*       std::cout<<"encountered nan in dt whil doing apply update"<<std::endl; */

/*     for(int i=0;i<6;i++) */
/*       { */
/*       if(isnan(  (u.Value())[i]   ))  */
/*         std::cout<<"encountered nan in delta value before dt multiply"<<std::endl; */
/*       if(isinf((u.Value())[i]) || isinf(-u.Value()[i])  )  */
/*         std::cout<<"encountered +/- inf in delta value before dt multiply"<<std::endl; */
/*       } */
    
    deltaValue= static_cast<PixelType>(u.Value() * dt);


    for(int i=0;i<6;i++)
      {
      refTensor[i]=originalValue[i];
      deltaTensor[i]=deltaValue[i];
      }

   if(IsValidTensor(refTensor)) 
    updatedTensor=ssTensorGeometry.ExpMap(refTensor,deltaTensor);
    bool isZero=true;
    bool flag=false;
	    
    for(int i=0;i<6;i++)
      {
      updatedValue[i]=updatedTensor[i];
      if(updatedValue[i]>10)
        flag=true;
      if(updatedValue[i]!=0)
        isZero=false;
      }

    if(isZero)
      count++;
      
    if(flag)
      {
    //  std::cout<<"refTensor="<<refTensor<<std::endl;
    //  std::cout<<"delta="<<deltaTensor<<std::endl;
    //  std::cout<<"updatedValue="<<updatedValue<<std::endl;
      }

    if(IsValidTensor(refTensor))
	    o.Value()=updatedValue ;  // no adaptor support here
    else
    {
	    invalid++;
    }
    ++o;
    ++u;
    
    }
    std::cout<<"After Threaded Apply update Zero Tensors="<<count<<std::endl;
    std::cout<<"Invalid count="<<invalid<<std::endl;
}

}//end namespace itk

#endif
