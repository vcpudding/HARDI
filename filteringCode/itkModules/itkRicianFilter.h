/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVectorGradientAnisotropicDiffusionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2003/09/10 14:28:59 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRicianFilter_h_
#define __itkRicianFilter_h_

#include "itkExceptionObject.h"
#include "itkAnisotropicDiffusionImageFilter.h"
#include "RicianAttachmentTerm.h"
#include "GaussianAttachmentTerm.h"
#include "itkDenseFiniteDifferenceImageFilter.h"
#include "itkVectorGradientNDAnisotropicDiffusionFunction.h"

namespace itk {

	
template <class TInputImage, class TOutputImage>
class ITK_EXPORT RicianFilter
  : public AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef RicianFilter Self;
  typedef AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>
  Superclass;
  typedef DenseFiniteDifferenceImageFilter<TInputImage, TOutputImage>
    DenseFiniteDifferenceFilterclass;

  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  
  /** The pixel type of the output image will be used in computations.
   * Inherited from the superclass. */
  typedef typename Superclass::PixelType PixelType;

  typedef typename PixelType::ValueType ScalarValueType;
  
  /** The value type of a time step. */
  typedef typename DenseFiniteDifferenceFilterclass::TimeStepType TimeStepType;
  
  /* The type of out put Image for this filter */
  typedef typename Superclass::OutputImageType OutputImageType;
  
  
  /** The container type for the update buffer. */
  typedef OutputImageType UpdateBufferType;
 
 
  /** The type of region used for multithreading */
  typedef typename UpdateBufferType::RegionType ThreadRegionType;

  typedef typename TInputImage::Pointer InputImagePointer;
   /** Instantiation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information. */
  itkTypeMacro(RicianFilter,
               AnisotropicDiffusionImageFilter);
  
  /** Determine the image dimension from the  superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension );

public:
  void SetNoisyImage(InputImagePointer noisyImage)
  {
	 this->noisyImage=noisyImage; 
  }
  
  void InitializeAttachmentTermObjects(double sigma,double lamda1,double lamda2)
  {
      rat= new RicianAttachmentTerm(sigma,lamda1);
      gat= new GaussianAttachmentTerm(sigma,lamda2);
      at=rat;//Default is a Rician Attachment Term
  }
protected:
  RicianFilter():at(0),rat(0),gat(0)
  {
    typename VectorGradientNDAnisotropicDiffusionFunction<UpdateBufferType>::Pointer p        
      = VectorGradientNDAnisotropicDiffusionFunction<UpdateBufferType>::New();
    this->SetDifferenceFunction(p);
  }
  ~RicianFilter() {
    delete gat;
    delete rat;
  }

private:
  typename TInputImage::Pointer noisyImage;  
  AttachmentTerm* at;
  RicianAttachmentTerm* rat;
  GaussianAttachmentTerm* gat;
  RicianFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  /**  Does the actual work of updating the output from the UpdateContainer over
   *  an output region supplied by the multithreading mechanism.
   *  Overloaded the ApplyUpdate in DenseFiniteDifferenceImageFilter
   *  to support update in Symmetric Space.
   */ 
  virtual
  void ThreadedApplyUpdate(TimeStepType dt,
                           const ThreadRegionType &regionToProcess,
                           int threadId);
};

} // end namspace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRicianFilter.txx"
#endif


#endif
