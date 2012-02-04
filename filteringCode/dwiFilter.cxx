#include<iostream>
#include<iomanip>
#include<sstream>

#include "itkNrrdImageIO.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkVariableLengthVector.h"         
#include "itkImageDuplicator.h"

#include "itkModules/itkRicianFilter.h"
#include "itkModules/itkGaussianFilter.h"

#define NUMIMAGES 7 
//#define NUMIMAGES 19 

typedef itk::Image<double, 3> DiffusionDataImageType;
typedef itk::ImageFileWriter<DiffusionDataImageType> DiffusionDataImageWriterType;
typedef itk::ImageFileReader<DiffusionDataImageType> DiffusionDataImageReaderType;
typedef itk::Vector< double,NUMIMAGES > VectorPixelType;
typedef itk::Image< VectorPixelType,3> VectorDiffusionDataImageType;
typedef itk::ImageDuplicator< DiffusionDataImageType> ImageDuplicator;
typedef itk::ImageDuplicator< VectorDiffusionDataImageType> ImageDuplicator2;

typedef itk::ImageFileWriter<VectorDiffusionDataImageType> VectorDiffusionDataImageWriterType;
typedef itk::ImageFileReader<VectorDiffusionDataImageType> VectorDiffusionDataImageReaderType;

typedef unsigned int uint32;

/*
 * createDWIVectorImage: returns a 3D volume of vectors whose components are
 * the dwi image signals along each gradient directions.NumImages specifies the
 * vector dimension.
 * dwi is an array of pointers to each 2D Diffusion Weighted Image
 * numImages: The number of 2D DWI images in the array.
 */
VectorDiffusionDataImageType::Pointer createDWIVectorImage(DiffusionDataImageType::Pointer dwi[],uint32 numImages)
{

	//Allocate Memory For the 3D Volume of vectors
	VectorDiffusionDataImageType::Pointer vectorDWI=VectorDiffusionDataImageType::New(); 
	VectorDiffusionDataImageType::IndexType start;
	start[0]=0;
	start[1]=0;
	start[2]=0;
	//define the size of the image
	DiffusionDataImageType::SizeType dwiSize= 
		dwi[0]->GetLargestPossibleRegion().GetSize();
	VectorDiffusionDataImageType::SizeType size;
	size[0]=dwiSize[0];
	size[1]=dwiSize[1];
	size[2]=dwiSize[2];
	//Initialize a region 
	VectorDiffusionDataImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);
	//Allocate memory
	vectorDWI->SetRegions(region);
	vectorDWI->SetSpacing(dwi[0]->GetSpacing());
	vectorDWI->SetDirection(dwi[0]->GetDirection());
	vectorDWI->Allocate();

	VectorDiffusionDataImageType::IndexType vectorPixelIndex;
	VectorDiffusionDataImageType::PixelType vectorPixelValue;
	DiffusionDataImageType::IndexType scalarPixelIndex;
	DiffusionDataImageType::PixelType scalarPixelValue;

	for(int i=0;i<size[0];i++) 
		for(int j=0;j<size[1];j++)
			for(int k=0;k<size[2];k++)
			{
				vectorPixelIndex[0]=i;
				vectorPixelIndex[1]=j;
				vectorPixelIndex[2]=k;

				scalarPixelIndex[0]=i;
				scalarPixelIndex[1]=j;
				scalarPixelIndex[2]=k;

				//Set the vector components
				for(unsigned int vdim=0;vdim<numImages;vdim++)
				{
					scalarPixelValue=dwi[vdim]->GetPixel(scalarPixelIndex);
					vectorPixelValue[vdim]=scalarPixelValue;
				}
				vectorDWI->SetPixel(vectorPixelIndex,vectorPixelValue);
			}
	return(vectorDWI);
}

/*
 * Read Real DWI data
 */
VectorDiffusionDataImageType::Pointer readDWI_ImagesREAL(char* prefix)
{
	VectorDiffusionDataImageReaderType::Pointer dwiReader= 
		VectorDiffusionDataImageReaderType::New();
	dwiReader->SetFileName(prefix);
	try{
		dwiReader->Update();
	}
	catch (itk::ExceptionObject e) 
	{
		std::cerr << e << std::endl;
		exit(-1);
	};
	return(dwiReader->GetOutput());
}


/*
 * Write out DWI data
 */
void writeDWI_ImagesREAL(char* prefix,
		VectorDiffusionDataImageType::Pointer image)
{
	
	VectorDiffusionDataImageWriterType::Pointer dwiWriter= 
		VectorDiffusionDataImageWriterType::New();
	dwiWriter->SetFileName(prefix);
	dwiWriter->SetInput(image);
	try{
		dwiWriter->Update();
	}
	catch (itk::ExceptionObject e) 
	{
		std::cerr << e << std::endl;
		exit(-1);
	}
}
/*
 * readDWI_Images returns a pointer to a 3D volume of vectors with dimension
 * determined by numImages parameter.
 */
VectorDiffusionDataImageType::Pointer readDWI_Images(const uint32 numImages,char* prefix)
{
	DiffusionDataImageReaderType::Pointer dwiReader= DiffusionDataImageReaderType::New();
	DiffusionDataImageType::Pointer  dwi[numImages];
	ImageDuplicator::Pointer duplicatorFilter=ImageDuplicator::New();
	std::ostringstream ostr;

	for(int i=0;i<numImages;i++)
	{
		ostr<<prefix<<std::setfill('0')<<std::setw(3)<<i<<".nrrd";
		std::cout<<"Reading File="<<ostr.str()<<std::endl;
		dwiReader->SetFileName((ostr.str()).c_str());
		ostr.str("");
		try
		{
			dwiReader->Update();
			duplicatorFilter->SetInputImage(dwiReader->GetOutput()); 
			duplicatorFilter->Update();
			dwi[i]=duplicatorFilter->GetOutput();
		}
		catch (itk::ExceptionObject e) 
		{
			std::cerr << e << std::endl;
			exit(-1);
		}

	}
	VectorDiffusionDataImageType::Pointer dwiVector=createDWIVectorImage(dwi,numImages);
	return(dwiVector);
}

/*
 * Square each value in the original
 * DWI volume
 */
VectorDiffusionDataImageType::Pointer
squareDWIVolume(VectorDiffusionDataImageType::Pointer vectorDWI)
{
      std::cout<<"Squaring Image....."<<std::endl;	
	VectorDiffusionDataImageType::SizeType size= 
		vectorDWI->GetLargestPossibleRegion().GetSize();
	
	VectorDiffusionDataImageType::IndexType vectorPixelIndex;
	VectorDiffusionDataImageType::PixelType vectorPixelValue;

	for(int i=0;i<size[0];i++) 
		for(int j=0;j<size[1];j++)
			for(int k=0;k<size[2];k++)
			{
				vectorPixelIndex[0]=i;
				vectorPixelIndex[1]=j;
				vectorPixelIndex[2]=k;
				vectorPixelValue=
					vectorDWI->GetPixel(vectorPixelIndex);
				for(int i=0;i<vectorPixelValue.Size();i++)
					vectorPixelValue[i]*=vectorPixelValue[i];
				vectorDWI->SetPixel(vectorPixelIndex,vectorPixelValue);
			}
	return(vectorDWI);
}

/*
 * Perform Rician Bias Correction
 */
VectorDiffusionDataImageType::Pointer
biasCorrect(VectorDiffusionDataImageType::Pointer vectorDWI,double sigma)
{
	std::cout<<"Performing Bias Correction..."<<std::endl;
	VectorDiffusionDataImageType::SizeType size= 
		vectorDWI->GetLargestPossibleRegion().GetSize();

	VectorDiffusionDataImageType::IndexType vectorPixelIndex;
	VectorDiffusionDataImageType::PixelType vectorPixelValue;

	for(int i=0;i<size[0];i++) 
		for(int j=0;j<size[1];j++)
			for(int k=0;k<size[2];k++)
			{
				vectorPixelIndex[0]=i;
				vectorPixelIndex[1]=j;
				vectorPixelIndex[2]=k;

				vectorPixelValue=
					vectorDWI->GetPixel(vectorPixelIndex);
				for(int i=0;i<vectorPixelValue.Size();i++){
					vectorPixelValue[i]-=(2*sigma*sigma);
					if(vectorPixelValue[i]<=0)
						vectorPixelValue[i]=0;
					else
						vectorPixelValue[i]=sqrt(vectorPixelValue[i]);
				}
				vectorDWI->SetPixel(vectorPixelIndex,vectorPixelValue);
			}
	return(vectorDWI);

}

/*
 * This routine does anisotropic Filtering on the 3-D volume of 
 * DWI Images
 */
VectorDiffusionDataImageType::Pointer anisoFilter(
		VectorDiffusionDataImageType::Pointer vectorDWI,
		uint32 numIterations,
		double conductance,
		double timeStep,
		uint32 filterType,
		double sigma,
		double lamda1=1.0,
		double lamda2=1.0)
{
	typedef itk::VectorGradientAnisotropicDiffusionImageFilter
		< VectorDiffusionDataImageType, VectorDiffusionDataImageType >  PeronaMalikFilterType;
	typedef itk::RicianFilter
		< VectorDiffusionDataImageType, VectorDiffusionDataImageType >  RicianFilterType;
	typedef itk::GaussianFilter
		< VectorDiffusionDataImageType, VectorDiffusionDataImageType >  GaussianFilterType;
	ImageDuplicator2::Pointer duplicatorFilter=ImageDuplicator2::New();
	duplicatorFilter->SetInputImage(vectorDWI); 
	duplicatorFilter->Update();
	VectorDiffusionDataImageType::Pointer noisyVectorDWI=duplicatorFilter->GetOutput();

	VectorDiffusionDataImageType::Pointer input;
	VectorDiffusionDataImageType::Pointer filteredVolume;

	std::cout<<"k="<<conductance;
	std::cout<<",dt="<<timeStep;
	std::cout<<",numiter="<<numIterations<<std::endl;

	PeronaMalikFilterType::Pointer filter1;
	RicianFilterType::Pointer filter2;
	GaussianFilterType::Pointer filter3;

	switch(filterType)
	{
		case 0:  //Simple Aniso Filtering
	
			std::cout<<"Simple AnisoTropic Filtering ..."<<std::endl;
			filter1= PeronaMalikFilterType::New();                		
			filter1->SetInput(vectorDWI);
			filter1->SetNumberOfIterations(numIterations);
			filter1->SetTimeStep(timeStep);
			filter1->SetConductanceParameter(conductance);
			filter1->Update();
			filteredVolume=(filter1->GetOutput()); 
			break;
		case 1: //Aniso Filtering squared image
			std::cout<<"Squared Magnitude Filtering.."<<std::endl;
			input=squareDWIVolume(vectorDWI);
			filter1=PeronaMalikFilterType::New();
			filter1->SetInput(vectorDWI);
			filter1->SetNumberOfIterations(numIterations);
			filter1->SetTimeStep(timeStep);
			filter1->SetConductanceParameter(conductance);
			filter1->Update();
			filteredVolume=biasCorrect(filter1->GetOutput(),sigma);	
			break;
		case 2: //Rician Correction Term
			std::cout<<"Rician Filtering"<<std::endl;
			std::cout<<"Bias Correction term"<<std::endl;
			filter2=RicianFilterType::New();
			filter2->SetInput(vectorDWI);
			filter2->SetNumberOfIterations(numIterations);
			filter2->SetTimeStep(timeStep);
			filter2->SetConductanceParameter(conductance);
			filter2->SetNoisyImage(vectorDWI);
			filter2->InitializeAttachmentTermObjects(sigma,lamda1,lamda2);
			filter2->Update();
			filteredVolume=(filter2->GetOutput()); 
			break;
		case 3: //Gaussian Correction Term
			std::cout<<"ANiso Filtering with Gaussian"<<std::endl;
			std::cout<<"Bias Correction term"<<std::endl;
			filter3=GaussianFilterType::New();
			filter3->SetInput(vectorDWI);
			filter3->SetNumberOfIterations(numIterations);
			filter3->SetTimeStep(timeStep);
			filter3->SetConductanceParameter(conductance);
			filter3->SetNoisyImage(noisyVectorDWI);
			filter3->InitializeAttachmentTermObject(sigma,lamda2);
			filter3->Update();
			filteredVolume=(filter3->GetOutput()); 
			break;

	}

	std::cout<<"Finished applying Filter"<<std::endl;

	return(filteredVolume);
}

/*
 * Convert vector DWI volume to scalar DWI volumes
 */
void createScalarDWIImages(VectorDiffusionDataImageType::Pointer vectorDWI,
		DiffusionDataImageType::Pointer dwi[],
		int numImages)
{

	//Allocate memory for numImage 2D DWI's
	DiffusionDataImageType::IndexType start;
	start[0]=0;
	start[1]=0;
	start[2]=0;
	//define the size of the image
	VectorDiffusionDataImageType::SizeType dwiSize= 
		vectorDWI->GetLargestPossibleRegion().GetSize();
	DiffusionDataImageType::SizeType size;
	size[0]=dwiSize[0];
	size[1]=dwiSize[1];
	size[2]=dwiSize[2];
	//Initialize a region 
	VectorDiffusionDataImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);

	//Allocate memory
	for(int i=0;i<numImages;i++)
	{
		dwi[i]=DiffusionDataImageType::New();
		dwi[i]->SetRegions(region);
		dwi[i]->SetSpacing(vectorDWI->GetSpacing());
		dwi[i]->SetDirection(vectorDWI->GetDirection());
		dwi[i]->Allocate();
	}

	VectorDiffusionDataImageType::IndexType vectorPixelIndex;
	VectorDiffusionDataImageType::PixelType vectorPixelValue;
	DiffusionDataImageType::IndexType scalarPixelIndex;
	DiffusionDataImageType::PixelType scalarPixelValue;

	for(int i=0;i<size[0];i++) 
		for(int j=0;j<size[1];j++)
			for(int k=0;k<size[2];k++)
			{
				vectorPixelIndex[0]=i;
				vectorPixelIndex[1]=j;
				vectorPixelIndex[2]=k;

				scalarPixelIndex[0]=i;
				scalarPixelIndex[1]=j;
				scalarPixelIndex[2]=k;

				//Set the vector components
				vectorPixelValue= vectorDWI->GetPixel(vectorPixelIndex);
				for(unsigned int vdim=0;vdim<numImages;vdim++)
					dwi[vdim]->SetPixel(scalarPixelIndex,vectorPixelValue[vdim]);
			}
}


/*
 * writeFilteredDWI writes out the vector DWI volume
 * as separate DWI Images which can then be estimated 
 * using the estimate tensor program
 */

void writeFilteredDWI(VectorDiffusionDataImageType::Pointer vectorDWI,
		const uint32 numImages,
		char *prefix)
{

	DiffusionDataImageType::Pointer dwi[numImages];
	createScalarDWIImages(vectorDWI,dwi,numImages);

	DiffusionDataImageWriterType::Pointer dwiWriter = DiffusionDataImageWriterType::New();
	std::ostringstream ostr;

	for(int i=0;i<numImages;i++)
	{
		ostr<<prefix<<std::setfill('0')<<std::setw(3)<<i<<".nrrd";
		std::cout<<"Writing File="<<ostr.str()<<std::endl;
		dwiWriter->SetFileName((ostr.str()).c_str());
		dwiWriter->SetInput(dwi[i]);
		ostr.str("");
		try
		{
			dwiWriter->Update();
		}
		catch (itk::ExceptionObject e) 
		{
			std::cerr << e << std::endl;
			exit(-1);
		}
	}
}


/*
 * Main Routine to perform
 * dwi filtering by squaring 
 * the dwi images
 */
int main(int argc, char* argv[])
{
	if(argc!=10)
	{
		std::cout<<"USAGE:dwiFilter <arguments>"<<std::endl;
		std::cout<<"Arguments:"<<std::endl;
		std::cout<<"1. Input File Name"<<std::endl;
		std::cout<<"2. Output File Name  "<<std::endl;
		std::cout<<"3. NumIterations"<<std::endl;
		std::cout<<"4. Conductance"<<std::endl;
		std::cout<<"5. TimeStep"<<std::endl;
		std::cout<<"6. Filter Type :";
		std::cout<<" (Simple Aniso-0,Chi Squared-1,Rician-2,Gaussian-3)"<<std::endl;
		std::cout<<"7. Sigma for bias correction"<<std::endl; 
		std::cout<<"8. Lamda (Rician Correction Term)"<<std::endl;
		std::cout<<"9. Lamda (Gaussian Correction Term)"<<std::endl; 
		exit(1);
	}
	else
	{
		uint32 numIter=static_cast<uint32>(atoi(argv[3])) ;
		double conductance=atof(argv[4]);
		double timeStep=atof(argv[5]);
		uint32 filterType= atoi(argv[6]);
		double sigma=atof(argv[7]);
		double lamda1=atof(argv[8]);
		double lamda2=atof(argv[9]);

		uint32 numImages=NUMIMAGES;
		VectorDiffusionDataImageType::Pointer vectorDWI=readDWI_ImagesREAL(argv[1]);

		VectorDiffusionDataImageType::Pointer smoothedvectorDWI=
			anisoFilter(vectorDWI,
				     numIter,
				     conductance,
				     timeStep,
				     filterType,
				     sigma,
				     lamda1,
				     lamda2);
			writeDWI_ImagesREAL(argv[2],smoothedvectorDWI);     
	}
	return (EXIT_SUCCESS);
}
