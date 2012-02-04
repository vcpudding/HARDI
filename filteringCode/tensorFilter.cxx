#include<iostream>
#include<iomanip>
#include "itkNrrdImageIO.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkDiffusionTensor3D.h"
#include "itkImageLinearIteratorWithIndex.h"
//Itk Module includes
#include "itkModules/itkSymmetricSpaceTensorAnisotropicDiffusionFilter.h"
#include "itkModules/itkEuclideanVectorToDiffusionTensorImageFilter.h"
#include "itkModules/itkDiffusionTensorToEuclideanVectorImageFilter.h"
#include "itkModules/itkGaussianFilter.h"
#include "itkImageDuplicator.h"
#include "itkModules/itkRicianFilter.h"

//Namic Includes-------------------
#include "SymmetricSpaceTensorGeometry.h"
#include "LinearTensorGeometry.h"
#include "TensorStatistics.h"

/*-----------------------------------------------------------
 * Driver Program for performing AnisoDiffusion 
 * on Diffusion Tensors  
 * Part of NAMIC Riemannian DTI Research Project 
 * Author:
 * Saurav Basu
 * SCI Institute
 * Start Date: 11/1/2005
 * Last Modified: 1/26/2006
 * University of Utah
 * NAMIC Project
 * --------------------------------------------------------*/

//ITK typedefs
typedef double PixelType;
typedef itk::CovariantVector< PixelType, 6 > VectorPixelType;
typedef itk::Image<VectorPixelType,3 > VectorImageType;
typedef itk::DiffusionTensor3D<PixelType> DiffusionTensorPixelType;
typedef itk::Image<DiffusionTensorPixelType, 3> DiffusionTensorImageType;
enum FilterSpaceType { EUCLIDEAN =0,LOGEUCLIDEAN,RIEMANNIAN};
typedef itk::Image<PixelType,3 > MaskImageType;
typedef itk::ImageDuplicator< VectorImageType> ImageDuplicator;

//Tensor Library typedefs
typedef itk::SymmetricSecondRankTensor<PixelType, 3> TensorType;


//helper method
template<class T>
bool IsValidTensor(T &tensor)
{
	//First check the reference tensor 
	bool valid=false;
	for(int i=0;i<6;i++)
		if(tensor[i])
			valid=true;
	return valid;
}
/*------------------------------------------
 * Check for zero/negative values in 
 * the TensorData 
 *-------------------------------------------*/
bool fixTensorData(VectorPixelType&  pixelValue,
		unsigned long* zeroCount,
		unsigned long* negativeCount
		)
{
	typedef DiffusionTensorPixelType::EigenValuesArrayType EigenValuesArrayType;
	typedef DiffusionTensorPixelType::EigenVectorsMatrixType EigenVectorsMatrixType;

	DiffusionTensorPixelType temp;

	for(int i=0;i<6;i++)
		temp[i]=pixelValue[i];

	bool fixed=false;
	bool isZero=true;
	for(unsigned int i=0;i<6;i++)
		if(temp[i]!=0)
			isZero=false;

	unsigned int indices[]={0,3,5};
	//double epsilon=1.0e-1;  //For torus data
	//double epsilon=1.0e-3; //For real data in paper
	//double epsilon =1.0e-4; //For synthetic (data paper)
	double epsilon=1.0e-5; //real data post proposal


	if(isZero)
	{
		//std::cout<<"pixelValue is Zero"<<std::endl;
		//for(unsigned int i=0;i<3;i++)
			//temp[indices[i]]=epsilon; 
		(*zeroCount)++;
		fixed=true;
	}
	else
	{
		EigenVectorsMatrixType U;
		EigenValuesArrayType  eVA;
		temp.ComputeEigenAnalysis(eVA,U);
	        double trace=temp.GetTrace(); 	
		double min=FLT_MAX;
		double max=FLT_MIN;
		for(int i=0;i<3;i++)
		{
			if(eVA[i]<=min)
				min=eVA[i];
			if(eVA[i]>=max)
				max=eVA[i];
		}
		if(min<=0 )
		{
			for(int i=0;i<3;i++)
				if(eVA[i]<0)
					eVA[i]+=(epsilon-min);
			vnl_matrix_fixed<DiffusionTensorPixelType::ComponentType,3,3> 
				D(0.0),vnl_U(U.GetVnlMatrix()),Result(0.0);
			D(0,0)=eVA[0];
			D(1,1)= eVA[1];
			D(2,2)=eVA[2];
			Result= U.GetTranspose()*D*vnl_U;
			temp[0]= Result(0,0);
			temp[1]= Result(0,1);
			temp[2]= Result(0,2);
			temp[3]= Result(1,1);
			temp[4]= Result(1,2);
			temp[5]= Result(2,2);
			fixed=true;
			(*negativeCount)++;
		}
	}
		for(int i=0;i<6;i++)
			pixelValue[i]=temp[i];
	
		return(fixed);
}
/*---------------------------------------------------
 *This function  fixes data problems
 * to avoid negative eigen values
 *---------------------------------------------------*/
VectorImageType::Pointer DataFix(VectorImageType::Pointer inputImage,
		int volWidth,
		int volHeight,
		int volDepth)
{
	
	VectorImageType::PixelType pixelValue,temp;
	VectorImageType::IndexType pixelIndex;
	
	bool fixed=false;
	unsigned long zeroCount=0,negativeCount=0; 
	
	
	for(int i=0;i<volWidth;i++) 
		for(int j=0;j<volHeight;j++)
			for(int k=0;k<volDepth;k++)
			{
				pixelIndex[0]=i;
				pixelIndex[1]=j;
				pixelIndex[2]=k;
				
				pixelValue=inputImage->GetPixel(pixelIndex);
				for(int m=0;m<6;m++)
					temp[m]=pixelValue[m];
				//Account for Negative Eigen Values 
				fixed=fixTensorData(temp,&zeroCount,&negativeCount);
				if(fixed)
					fixed=false;
				inputImage->SetPixel(pixelIndex,temp);
			}
	std::cout<<"Zero Count="<<zeroCount<<std::endl;
	std::cout<<"Negative Count="<<negativeCount<<std::endl;

	return(inputImage);
}


/*---------------------------------------------------
 *This function applies the log map 
 * to a vector image of tensors  
 *---------------------------------------------------*/
VectorImageType::Pointer applyMap(
		VectorImageType::Pointer inputImage,
		bool mapType, 
		int volWidth,
		int volHeight,
		int volDepth)
{
	TensorType inputTensor,outputTensor,mean;
	mean.SetIdentity();
	VectorImageType::PixelType pixelValue;
	VectorImageType::IndexType pixelIndex;

	//Allocate memory for output image
	VectorImageType::Pointer outputImage = VectorImageType::New();
	VectorImageType::IndexType start;
	start[0]=0;
	start[1]=0;
	start[2]=0;
	VectorImageType::SizeType size;
	size[0]=volWidth;
	size[1]=volHeight;
	size[2]=volDepth;
	VectorImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);
	outputImage->SetRegions(region);
	outputImage->SetSpacing(inputImage->GetSpacing());
	outputImage->SetDirection(inputImage->GetDirection());
	outputImage->Allocate();

	//Create geometry object
	SymmetricSpaceTensorGeometry<VectorPixelType::ValueType,3> * ss_Geom;
	ss_Geom = new SymmetricSpaceTensorGeometry<VectorPixelType::ValueType,3>;

	for(int i=0;i<volWidth;i++) 
		for(int j=0;j<volHeight;j++)
			for(int k=0;k<volDepth;k++)
			{
				pixelIndex[0]=i;
				pixelIndex[1]=j;
				pixelIndex[2]=k;
				pixelValue=inputImage->GetPixel(pixelIndex);
				for(int l=0;l<6;l++)
					inputTensor[l]=pixelValue[l];

				if(IsValidTensor(inputTensor))
				outputTensor=(mapType)?(ss_Geom->LogMap(mean,inputTensor))
					:(ss_Geom->ExpMap(mean,inputTensor));

				for(int l=0;l<6;l++)
					pixelValue[l]=outputTensor[l];
				outputImage->SetPixel(pixelIndex,pixelValue);
			}
	delete ss_Geom;
	return(outputImage);
}
/*---------------------------------------------------
 *This function applies the Log map 
 * to a vector image of tensors
 *---------------------------------------------------*/
VectorImageType::Pointer applyLogMap(VectorImageType::Pointer inputImage,
		int volWidth,
		int volHeight,
		int volDepth)
{
	return(applyMap(inputImage,true,volWidth,volHeight,volDepth));
}

/*---------------------------------------------------
 *This function applies the Exp map 
 * to a vector image of tensors
 *---------------------------------------------------*/
VectorImageType::Pointer applyExpMap(VectorImageType::Pointer inputImage,
		int volWidth,
		int volHeight,
		int volDepth)
{
	return(applyMap(inputImage,false,volWidth,volHeight,volDepth));
}
/*---------------------------------------------------
 * This function applies anisotropic diffusion 
 *  for vector valued images 
 *  FilterSpace= EUCLIDEAN/RIEMANNIAN/LOG SPACE
 *- -------------------------------------------------*/
VectorImageType::Pointer anisoDiffuse(VectorImageType::Pointer  inputImage,
		FilterSpaceType filterSpace,
		int volWidth,
		int volHeight,
		int volDepth,
		int numIterations=5,
		double timeStep=0.01,
		double conductance=1.0
		)
{

	std::cout<<"Num Iterations="<<numIterations<<std::endl;
	std::cout<<"Starting Anisotropic Diffusion on Volume of Size=( "
		<<volWidth<<","<<volHeight<<","<<volDepth<<")"<<std::endl;
	typedef itk::VectorGradientAnisotropicDiffusionImageFilter
		< VectorImageType, VectorImageType > EuclideanFilterType;
	typedef itk::GaussianFilter
		< VectorImageType, VectorImageType > GaussianFilterType;
	typedef itk::RicianFilter
		< VectorImageType, VectorImageType > RicianFilterType;

	typedef itk::SymmetricSpaceTensorAnisotropicDiffusionFilter
		< VectorImageType, VectorImageType > SymmetricSpaceFilterType;
	typedef itk::AnisotropicDiffusionImageFilter
		<VectorImageType ,VectorImageType  > FilterType;

	FilterType::Pointer filter;
	VectorImageType::Pointer filterInput;
	VectorImageType::Pointer filterOutput,temp;

	switch(filterSpace)
	{
		case EUCLIDEAN:
			filter=EuclideanFilterType::New();
			//Process the image using the filter
			filter->SetInput(inputImage);
			filter->SetNumberOfIterations(numIterations);
			filter->SetTimeStep(timeStep);
			filter->SetConductanceParameter(conductance);
			std::cout<<"Starting to apply Euclidean Space Filter with Gaussian Correction term..."<<std::endl;
			filter->Update();
			std::cout<<"Finished applying Euclidean Space Filter"<<std::endl;
			filterOutput=filter->GetOutput();
			break;
		case LOGEUCLIDEAN:
			filter=EuclideanFilterType::New();
			std::cout<<"Starting to  apply Log Map..."<<std::endl;
			filterInput=applyLogMap(inputImage,volWidth,volHeight,volDepth);
			std::cout<<"Finished applying Log Map"<<std::endl;
			filter->SetInput(filterInput);
			filter->SetNumberOfIterations(numIterations);
			filter->SetTimeStep(timeStep);
			filter->SetConductanceParameter(conductance);
			std::cout<<"Starting to apply Log Space Filter..."<<std::endl;
			filter->Update();
			std::cout<<"Finished applying Log Space Filter"<<std::endl;
			filterOutput=applyExpMap(filter->GetOutput(),volWidth,volHeight,volDepth);
			break;
		case RIEMANNIAN:
			filter=SymmetricSpaceFilterType::New();
			filter->SetInput(inputImage);
			filter->SetNumberOfIterations(numIterations);
			filter->SetTimeStep(timeStep);
			filter->SetConductanceParameter(conductance);
			filter->SetNumberOfThreads(1 );
			std::cout<<"Starting to apply Riemannian Space Filter..."<<std::endl;
			filter->Update();
			std::cout<<"Finished applying Riemannian Space Filter"<<std::endl;
			filterOutput=filter->GetOutput();
			break;
		default:
			std::cout<<"Invalid filter space:Terminating Program"<<std::endl;
			exit(1);
			break;
	}
	return(filterOutput);
}

/*-------------------------------------------------------------------
 * Read a Nrrd File of Diffusion Tensor Data
 *-------------------------------------------------------------------*/
DiffusionTensorImageType::Pointer readTensorNrrdFile( 
		char* fileName,
		int* volWidth,
		int* volHeight, 
		int* volDepth)
{
	itk::ImageFileReader<DiffusionTensorImageType>::Pointer reader
		= itk::ImageFileReader<DiffusionTensorImageType>::New();
	reader->SetFileName(fileName);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file reader " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		exit(-1);
	}
	DiffusionTensorImageType::Pointer image=(reader->GetOutput());
	DiffusionTensorImageType::SizeType imageSize= 
		image->GetLargestPossibleRegion().GetSize();

	(*volWidth)=imageSize[0];
	(*volHeight)=imageSize[1];
	(*volDepth)=imageSize[2];
	
	return(image);
}
/*------------------------------------------------------
 *This functions writes a VectorImage (of Tensors)  as  
 * a nrrd File
 *----------------------------------------------------*/
void writeTensorNrrdFile(char* fileName,
		DiffusionTensorImageType::Pointer image,
		int volWidth,
		int volHeight, 
		int volDepth)
{

	DiffusionTensorImageType::SizeType imageSize= 
		image->GetLargestPossibleRegion().GetSize();
	std::cout<<"Trying to write image of size ="<<imageSize<<std::endl;
	std::cout<<"Filename="<<fileName<<std::endl;
	itk::ImageFileWriter<DiffusionTensorImageType>::Pointer writer;
	writer = itk::ImageFileWriter<DiffusionTensorImageType>::New();
	writer->SetInput( image );
	writer->SetFileName(fileName);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "exception in file writer " << std::endl;
		std::cerr << e.GetDescription() << std::endl;
		std::cerr << e.GetLocation() << std::endl;
		exit(-1);
	}
}

/*----------------------------------------------------------
 *  Convert a Vector Image to a DiffusionTensor3D Image
 *--------------------------------------------------------*/
DiffusionTensorImageType::Pointer EuclideanVectorToDiffusionTensorImage(
		VectorImageType::Pointer & inputImage)
{
	typedef itk::EuclideanVectorToDiffusionTensorImageFilter<VectorImageType,
		DiffusionTensorImageType > WeightFilterType;
	
	WeightFilterType::Pointer weightFilter= WeightFilterType::New();
	weightFilter->SetInput(inputImage);
	weightFilter->Update();
	return(weightFilter->GetOutput());
}
/*-------------------------------------------------------------
 * Convert a  DiffusionTensor3D Image to VectorImage
 *-------------------------------------------------------------*/
VectorImageType::Pointer DiffusionTensorToEuclideanVectorImage(
		DiffusionTensorImageType::Pointer & inputImage)
{
	typedef itk::DiffusionTensorToEuclideanVectorImageFilter<DiffusionTensorImageType,
		VectorImageType> InverseFilterType;

	InverseFilterType::Pointer invFilter= InverseFilterType::New();
	invFilter->SetInput(inputImage);
	invFilter->Update();
	return(invFilter->GetOutput());
}

/*------------------------------------------------
 * Main Routine (Driver For the Filtering process)
 *--------------------------------------------------*/
int main(int argc , char* argv[])
{
	// Parameters
	double stdDev=0.2;
	int volWidth=10,volHeight=10,volDepth=10;
	int numIterations=30;
	double timeStep=0.01,conductance=1.0;
	bool useEuclideanSpaceNoise=false;
        double eigenValueAverage=0;
	FilterSpaceType filterSpace;
        double lamda=1.0;
	char filterNames[][20]={
		"Euclidean",
		"Log Space",
		"Riemannian"
	};

	char smoothed[][100]={
		"EuclideanFiltered.nhdr",
		"LogSpaceFiltered.nhdr",
		"SymmetricSpaceFiltered.nhdr"
	};


	if(argc!=6)
	{
		std::cout<<"Usage: tensorDiffuseTest <Arguments>"<<std::endl;
		std::cout<<"1. FilterType:(0-Euclidean, 1-Log Space,2-Riemannian)"<<std::endl;
		std::cout<<"2. numIterations:Iterations For Anisotropic Diffusion"<<std::endl;
		std::cout<<"3. timeStep:timeStep Used in Anisotropic Diffusion"<<std::endl;
		std::cout<<"4. conductance:Conductance used for Anisotropic Diffusion"<<std::endl;
		std::cout<<"5  Input (filename of input data)"<<std::endl;
		exit(1);
	}
	else
	{
		filterSpace=static_cast<FilterSpaceType>(atoi(argv[1]));
		numIterations=atoi(argv[2]);
		timeStep=atof(argv[3]);
		conductance=atof(argv[4]);

               if((filterSpace!=EUCLIDEAN) && (filterSpace!=RIEMANNIAN)
			       && (filterSpace!=LOGEUCLIDEAN))
	       {
		       std::cout<<"Invalid Filter space..Terminating"<<std::endl;
		       exit(1);
	       }
		
		std::cout<<" Performing Tensor "<<filterNames[filterSpace]
			<<" Filtering with parameters"<<std::endl;
		std::cout<<"Num Iterations="<<numIterations<<std::endl;
		std::cout<<"timeStep="<<timeStep<<std::endl;
		std::cout<<"conductance="<<conductance<<std::endl;
	}
	
	VectorImageType::Pointer noisyTensor,originalTensor,
				 smoothedTensor;
        
	DiffusionTensorImageType::Pointer originalNrrd,
				 noisyNrrd,
				 smoothedNrrd;
	MaskImageType::Pointer maskImage=NULL;
		//(Create/Read) Corrupted Tensor by adding  Noise to original Tensor volume 	
        	noisyNrrd=readTensorNrrdFile(argv[5],&volWidth,&volHeight,&volDepth); 
		noisyTensor=DiffusionTensorToEuclideanVectorImage(noisyNrrd);

		std::cout<<"Image size=("<<volWidth<<","<<volHeight<<","<<volDepth<<")"<<std::endl;
		noisyNrrd=EuclideanVectorToDiffusionTensorImage(noisyTensor); 

		//Log Space and Riemannian Filters require all tensors to be symmetric and 
		//positive definite. DATAFIX fixes the tensors which have negative eigen values
		// Euclidean space anisotropic filtering does not require this. 
		if((filterSpace!=EUCLIDEAN))
			noisyTensor=DataFix(noisyTensor,
					volWidth,
					volHeight,
					volDepth
					);   


        //Perform Filtering
	std::cout<<"Starting "<<filterNames[filterSpace]<< " Anisotropic Diffusion..."<<std::endl;
	smoothedTensor=anisoDiffuse(noisyTensor,
			filterSpace,
			volWidth,
			volHeight,
			volDepth,
			numIterations,
			timeStep,
			conductance);
	std::cout<<"Finished "<<filterNames[filterSpace]<<" Anisotropic Diffusion"<<std::endl;
	//Convert Vector Space back to tensors and write out results 
	smoothedNrrd=EuclideanVectorToDiffusionTensorImage(smoothedTensor);
	writeTensorNrrdFile(smoothed[filterSpace],smoothedNrrd,volWidth,volHeight,volDepth);
	return 0;
}

