#include "Riostream.h" 
#include "../interface/RooHighmass_conv.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"

using namespace TMath;

ClassImp(RooHighmass_conv) 

  RooHighmass_conv::RooHighmass_conv(const char *name, const char *title, 
					     RooAbsReal& _xreco,
//					     RooAbsReal& _xgen,
					     const RooArgList& inCoefList): 
   RooAbsPdf(name,title), 
   xreco("xreco","xreco",this,_xreco),
   //xgen("xgen","xgen",this,_xgen),
  _coefList("coefList","List of funcficients",this) 
  
 { 
  TIterator* coefIter = inCoefList.createIterator() ;
  RooAbsArg* func;
  while((func = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(func)) {
      coutE(InputArguments) << "ERROR: :RooHighmass_conv(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
      assert(0);
    }
    _coefList.add(*func) ;
  }
  delete coefIter;
	 RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(xreco.arg());
	 TIterator *arglist= xlist->createIterator();
	 RooRealVar *xgencopy= (RooRealVar*)arglist->Next();

	 xgen_min = xgencopy->getMin();
	 xgen_max = xgencopy->getMax();
  
 } 


 RooHighmass_conv::RooHighmass_conv(const RooHighmass_conv& other, const char* name) :  
   RooAbsPdf(other,name), 
   xreco("xreco",this,other.xreco),
//   xgen("xgen",this,other.xgen),
  _coefList("coefList",this,other._coefList)

 { 
	xgen_min = other.xgen_min;
	xgen_max = other.xgen_max;
 } 


 Double_t RooHighmass_conv::evaluate() const 
 { 
   double value = 0.;
//	 xgencopy->Print();
	 RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(xreco.arg());
	 TIterator *arglist= xlist->createIterator();
	 RooRealVar *xgencopy= (RooRealVar*)arglist->Next();

	for(int i=0;i<1000;i++){
		//dynamic_cast<RooRealVar >(xgen.arg()).setVal(i);
//		setVal(xgen_min+ i* (xgen_max-xgen_min)/100.);
	//xgen(xgen_min);
	 xgencopy->setVal(xgen_min+  i*(xgen_max-xgen_min)/1000.);
   Double_t reso = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(xreco.arg());
   Double_t gen = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
	 //cout << xgen_min+  i*(xgen_max-xgen_min)/1000.<<" "<<gen<<" "<< reso<<endl;
	 value += reso*gen;
	}
	return value;
   
 } 

Int_t RooHighmass_conv::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,xreco)) return 4 ;
//	else return 3;
//  if (matchArgs(allVars,analVars,xgen)) return 3 ;

  return 0 ;

}

Double_t RooHighmass_conv::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

     case 4: 
       {
		double integral=0;
	 RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(xreco.arg());
	 TIterator *arglist= xlist->createIterator();
	 RooRealVar *xgencopy= (RooRealVar*)arglist->Next();
	 //double xgen_min = xgencopy->getMin();
	 //double xgen_max = xgencopy->getMax();
		for(int j=0;j<1000;j++){
		 xgencopy->setVal(xgen_min+  j*(xgen_max-xgen_min)/1000.);
//  	 Double_t reso = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getNorm(xreco.arg());
   		Double_t gen = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
	 //		cout << gen<<" "<< reso<<endl;
	 	integral += gen;
		}
//		cout << integral<<endl;
		 return integral;

       }
       
     }
   
   assert(0) ;
   return 0 ;
}

