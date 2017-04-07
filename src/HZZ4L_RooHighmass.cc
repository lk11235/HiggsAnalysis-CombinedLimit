#include "Riostream.h" 
#include "../interface/HZZ4L_RooHighmass.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"

using namespace TMath;
using namespace std;

ClassImp(HZZ4L_RooHighmass) 

  HZZ4L_RooHighmass::HZZ4L_RooHighmass(const char *name, const char *title, 
					     RooAbsReal& _mass,
					     RooAbsReal& _mean,
					     RooAbsReal& _width,
					     RooAbsReal& _coupl,
					     const RooArgList& inCoefList): 
   RooAbsPdf(name,title), 
   mass("mass","mass",this,_mass),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   coupl("coupl","coupl",this,_coupl),
 // _normSet("normSet","List of observables",this), 
  _coefList("coefList","List of funcficients",this) 
  
 { 
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :HZZ4L_RooHighmass(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;
  
 } 


 HZZ4L_RooHighmass::HZZ4L_RooHighmass(const HZZ4L_RooHighmass& other, const char* name) :  
   RooAbsPdf(other,name), 
//	 _normSet(other._normSet),
   mass("mass",this,other.mass),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   coupl("coupl",this,other.coupl),
  //_normSet("normSet",this,other._normSet),
  _coefList("coefList",this,other._coefList)

 { 
  _coefIter = _coefList.createIterator() ;
 } 


 Double_t HZZ4L_RooHighmass::evaluate() const 
 { 
   double value = 0.;

//	 std::cout <<"xiaomeng"<<std::endl;	
//   Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal();
   Double_t T1 = dynamic_cast<const RooHistFunc*>(_coefList.at(0))->getVal();
   Double_t T2 = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->getVal();
	 //
 //  Double_t T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();
 
   RooArgSet *nset = new RooArgSet(mass.arg());
//   Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(nset);
//   Double_t T2 = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(nset);
   //Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal(nset)/5.;
   Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal(nset);
//	 cout << mass<<" "<<T3<<endl;
 //  Double_t T3 = dynamic_cast<const RooHistFunc*>(_coefList.at(2))->getVal();
//   Double_t T4 = dynamic_cast<const RooHistFunc*>(_coefList.at(3))->getVal();

//	 double a= (mass*mass-mean*mean);
//	 double b = mean*width;
//   double cossig = a/TMath::Sqrt(a*a+b*b);
//   double sinsig = b/TMath::Sqrt(a*a+b*b);


//	 std::cout <<"xiaomeng "<< T1<<" "<<T2<<std::endl;
//	 std::cout <<"xiaomeng "<<mass.arg().getVal()<<" "<< T2<<std::endl;

   
//   value = nbkg*T2 + fabs(coupl) * nsig*T1 + 1.76*sqrt(T2*T1*nbkg*nsig*coupl)*(cossig*T3+sinsig*T4); 
   //value = nbkg*T2 + fabs(coupl) * nsig*T1 + 1.76*sqrt(nbkg*nsig)*T3;
   value = T2 + coupl*T1 + sqrt(coupl)*T3;
	 //cout << mass<<" "<< value<<endl;
   
   if ( value <= 0.) { return 1.0e-200;}
   
   return value ; 
   
 } 

Int_t HZZ4L_RooHighmass::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,mass)) return 4 ;

  return 0 ;

}

Double_t HZZ4L_RooHighmass::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

     case 4: 
       {
//				 RooArgSet _norm
// double Int_T1  = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))-> getNorm(nset);
// double Int_T2  = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))-> getNorm(nset);
// double Int_T1 = dynamic_cast<const RooExtendPdf*>(_coefList.at(0))->expectedEvents(nset); 
// double Int_T2 = dynamic_cast<const RooExtendPdf*>(_coefList.at(1))->expectedEvents(nset); 
 double Int_T1  = dynamic_cast<const RooHistFunc*>(_coefList.at(0))-> analyticalIntegral(1000);
 double Int_T2  = dynamic_cast<const RooHistFunc*>(_coefList.at(1))-> analyticalIntegral(1000);
 //cout << "xiaomeng "<< dynamic_cast<const RooExtendPdf*>(_coefList.at(1))->expectedEvents(nset)<< " "<<  dynamic_cast<const RooExtendPdf*>(_coefList.at(0))->expectedEvents(nset)<<endl ;
// double Int_T2  = dynamic_cast<const RooHistFunc*>(_coefList.at(1))->dataHist().sum(false); 
// cout << nsig<< " "<< nbkg<<endl;
// double Int_T4  = dynamic_cast<const RooHistFunc*>(_coefList.at(2))-> analyticalIntegral(1000);

//   Double_t T3_int= dynamic_cast<const RooHistFunc*>(_coefList.at(2))->dataHist().sum(false);
   Double_t T3_int= dynamic_cast<const RooHistFunc*>(_coefList.at(2))->analyticalIntegral(1000);


	 double mysgn = 1.;
	 if(coupl < 0.) 
	   {
	     mysgn = -1.;
	   }

	 double integral =  coupl*Int_T1 + Int_T2 + mysgn*sqrt(fabs(coupl)) * T3_int; ;
//	 double integral =   Int_T1 + fabs(coupl) * Int_T2 + mysgn*sqrt(fabs(coupl)) * Int_T2 * 2 *(pow(mass,2)-pow(mean,2)) * Int_T1; 
//	    value = nbkg*T2 + fabs(coupl) * nsig*T1 + 1.76*sqrt(T2*T1*nbkg*nsig*coupl)*(cossig*T3+sinsig*T4);
//

	 //double integral =   nbkg + fabs(coupl) * nsig + sqrt(nbkg*nsig)*1.76*T3_int ;
	 if ( integral<= 0.) { return 1.0e-200;}
	 return integral;
       }
       
     }
   
   assert(0) ;
   return 0 ;
}

