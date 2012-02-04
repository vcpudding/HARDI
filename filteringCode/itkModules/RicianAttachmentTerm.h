#ifndef __RICIANATTACHMENT_TERM_H
#define __RICIANATTACHMENT_TERM_H
#include "AttachmentTerm.h"
#include "../besselFunction.h"

class RicianAttachmentTerm : public AttachmentTerm
{ 
public:	
   RicianAttachmentTerm(double sigma,double lamda):AttachmentTerm(sigma),_lamda(lamda){};
   double GetAttachmentTerm(double u0,double u)
   {
	   double bias=0;
	   double sigmaSquared=(_sigma*_sigma);	   
	   double term=u0/(sigmaSquared);
           double besselCutoff=700;
	   
       if(( term*u) < besselCutoff )
       {
          bias= (term*bessi1(term*u)/bessi0(term*u))- u/sigmaSquared;	       	
       //std::cout<<"Bias="<<bias<<std::endl;
       //std::cout<<"term="<<term<<",u="<<u<<",u0="<<u0<<"sigma="<<_sigma<<std::endl;       
       }
       else
       {
	       bias= term-u/sigmaSquared;
       }
	return (bias)*_lamda;
   }
private:
   double _lamda;
};

#endif
