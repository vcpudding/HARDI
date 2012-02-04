#ifndef __ATTACHMENT_TERM_H
#define __ATTACHMENT_TERM_H

class AttachmentTerm
{
public:	
   AttachmentTerm(double sigma):_sigma(sigma){};
   ~AttachmentTerm(){};	
   /*
    * u0 is noisy image
    * u is image being filtered
    */
   virtual double GetAttachmentTerm(double u0,double u){ return 0.0;}
   double Getsigma(){ return(_sigma);}
   
protected:
   double _sigma;
};
#endif
