#ifndef __GAUSSIANATTACHMENT_TERM_H
#define __GAUSSIANATTACHMENT_TERM_H

#include "AttachmentTerm.h"

class GaussianAttachmentTerm : public AttachmentTerm
{ 
public:	
   GaussianAttachmentTerm(double sigma,double lamda):AttachmentTerm(sigma),_lamda(lamda){};
   double GetAttachmentTerm(double u0,double u)
   {
        double bias=0;
       double sigmaSquared=(_sigma*_sigma);	   
       bias= 1/sigmaSquared*(u0-u);	       	
       //std::cout<<"Bias="<<bias<<std::endl;
       //std::cout<<"u="<<u<<",u0="<<u0<<"sigma="<<_sigma<<std::endl;       
       //std::cout<<"(u,u0)=("<<u<<","<<u0<<")"<<std::endl;
       //std::cout<<"lamda="<<_lamda<<std::endl; 
       //std::cout<<"bias="<<bias<<std::endl;
       return ((bias)*_lamda);
   }
private:
   double _lamda;
};

#endif
