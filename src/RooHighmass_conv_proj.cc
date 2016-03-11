#include "Riostream.h" 
#include "../interface/RooHighmass_conv_proj.h" 
#include "RooAbsReal.h" 
#include "RooProjectedPdf.h" 
#include "RooProdPdf.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"

using namespace TMath;

ClassImp(RooHighmass_conv_proj) 

  RooHighmass_conv_proj::RooHighmass_conv_proj(const char *name, const char *title, 
	 RooProdPdf &_intpdf, RooRealVar &_intObs,const RooArgSet &mreco
					     ): 
   RooProjectedPdf(name,title,_intpdf,_intObs), 
  intpdf("intpdf","pdfs",this,_intpdf), 
  _pdfList("pdfList","List of pdfs",this), 
  intObs("intobs","obs",this,_intObs), 
  xlist("xlist","xlist",this)
//  _coefList("coefList","List of funcficients",this) 
  
 { 
//  TIterator* coefIter = intObs.createIterator() ;
//  RooAbsArg* func;
//  while((func = (RooAbsArg*)coefIter->Next())) {
//    if (!dynamic_cast<RooAbsReal*>(func)) {
//      coutE(InputArguments) << "ERROR: :RooHighmass_conv_proj(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
//      assert(0);
//    }
//    _coefList.add(*func) ;
//		xgen_min = dynamic_cast<RooRealVar*>(func)->getMin();
//		xgen_max = dynamic_cast<RooRealVar*>(func)->getMax();
//  }
//  delete coefIter;
	//coefIter = _ciefList.createIterator();
	xgen_min = _intObs .getMin();
	xgen_max = _intObs .getMax();
//	_intObs.Print();
//	cout<< xgen_min<<" "<<xgen_max<<endl;

	RooArgList pdflist=_intpdf.pdfList();
	_pdfList.add(pdflist);

	 //xlist.add(*(_intpdf.getParameters(_coefList)));
	 xlist.add(mreco);
	 xlist.Print();

	 //TIterator *arglist= xlist->createIterator();
//	 xreco= (RooRealVar*)arglist->Next();
 } 


 RooHighmass_conv_proj::RooHighmass_conv_proj(const RooHighmass_conv_proj& other, const char* name) :  
   RooProjectedPdf(other,name), 
  intpdf("intpdf",this,other.intpdf),
  _pdfList("pdfList",this,other._pdfList),
  //_coefList("coefList",this,other._coefList)
  intObs("intobs",this,other.intObs),
  xlist("xlist",this,other.xlist)

 { 
	 xgen_min = other.xgen_min;
	 xgen_max = other.xgen_max;
 // _coefIter = _coefList.createIterator() ;
 } 


 Double_t RooHighmass_conv_proj::evaluate() const 
 { 
   double value = 0.;

   RooAbsPdf *pdf1 = (RooAbsPdf*)_pdfList.at(0);
   RooAbsPdf *pdf2 = (RooAbsPdf*)_pdfList.at(1);
//	 cout <<"xiaomeng"<<endl;
//	 _coefList.Print();
//	 RooAbsReal *integral =intpdf.arg().createIntegral(_coefList);
	 //cout << integral->getVal()<<endl;
//	 return integral->getVal();

	 RooArgSet *tmplist = pdf1->getParameters(xlist);
   TIterator *arglist= tmplist->createIterator();
   RooRealVar *mzz= (RooRealVar*)arglist->Next();

	for(int i=0;i<1000;i++){
	 double val = xgen_min+  i*(xgen_max-xgen_min)/1000.;
	 mzz->setVal(val);
	// cout << xgen_max <<" "<< xgen_min<<" " <<mzz->getVal()<<endl;
   Double_t reso = pdf2->getVal(xlist);
//	 cout << reso<<endl;
   Double_t gen = pdf1->getVal(*mzz);
//	 cout << reso<<" "<< gen<<endl;
	 value += reso*gen;
	}
	return value;
   
 } 

Int_t RooHighmass_conv_proj::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if (matchArgs(allVars,analVars,xlist)) return 4 ;

  return 0 ;

}

Double_t RooHighmass_conv_proj::analyticalIntegral(Int_t code, const char* rangeName) const
{
   switch(code)
     {

     case 4: 
       {
		double integral=0;

   RooAbsPdf *pdf1 = (RooAbsPdf*)_pdfList.at(0);
//   RooAbsPdf *pdf2 = (RooAbsPdf*)_pdfList.at(1);

	 RooArgSet *tmplist = pdf1->getParameters(xlist);
   TIterator *arglist= tmplist->createIterator();
   RooRealVar *mzz= (RooRealVar*)arglist->Next();

//	 double xgen_min = _coefList.at(0)->getMin();
//	 double xgen_max = _coefList.at(0)->getMax();

		for(int j=0;j<1000;j++){
		 double val = xgen_min+  j*(xgen_max-xgen_min)/1000.;
 		 mzz->setVal(val);
   		Double_t gen = pdf1->getVal();
	 	integral += gen;
		}
		 return integral;

       }
       
     }
   
   assert(0) ;
   return 0 ;
}

