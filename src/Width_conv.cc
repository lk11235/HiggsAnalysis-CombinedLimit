#include "Riostream.h" 
#include "RooDataSet.h" 
#include "../interface/Width_conv.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h>
#include "TMath.h"
#include "TH3F.h"
#include "TAxis.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TF1.h"
#include "TFile.h"
#include "TSpline.h"

using namespace TMath;

ClassImp(Width_conv) 

	Width_conv::Width_conv(const char *name, const char *title, 
			RooAbsReal& _xreco,
			RooAbsReal& _mean,
			RooAbsReal& _width,
			RooAbsReal& _coupl,
			const RooArgList& inCoefList, TGraph _cosspline, TGraph _sinspline, TGraph effsig, TGraph effbkg): 
		RooAbsPdf(name,title), 
		xreco("xreco","xreco",this,_xreco),
		mean("mean","mean",this,_mean),
		width("width","width",this,_width),
		coupl("coupl","coupl",this,_coupl),
		_coefList("coefList","List of funcficients",this) 

{ 
	TIterator* coefIter = inCoefList.createIterator() ;
	RooAbsArg* func;
	while((func = (RooAbsArg*)coefIter->Next())) {
		if (!dynamic_cast<RooAbsReal*>(func)) {
			coutE(InputArguments) << "ERROR: :Width_conv(" << GetName() << ") funcficient " << func->GetName() << " is not of type RooAbsReal" << endl;
			assert(0);
		}
		_coefList.add(*func) ;
	}
	delete coefIter;
	RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(mean.arg());
	TIterator *arglist= xlist->createIterator();
	xlist->Print();
	RooRealVar* xgencopy= (RooRealVar*)arglist->Next();

	RooArgSet *ylist= dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getParameters(*xgencopy);
	TIterator *arglisty= ylist->createIterator();
	RooRealVar* xrecocopy= (RooRealVar*)arglisty->Next();

	xgen_min = xgencopy->getMin();
	xgen_max = xgencopy->getMax();

	xreco_min = xrecocopy->getMin();
	xreco_max = xrecocopy->getMax();

	cout << xgen_max <<"\t"<<xgen_min<<endl;
	cout << xreco_max <<"\t"<<xreco_min<<endl;

	nsig= dynamic_cast<const RooExtendPdf*>(_coefList.at(0))->expectedEvents(*xgencopy); 
	nbkg= dynamic_cast<const RooExtendPdf*>(_coefList.at(1))->expectedEvents(*xgencopy); 
	//nbkg= 0; 
	cosspline=_cosspline;
	sinspline=_sinspline;
	cout << nsig<<"\t"<<nbkg<<endl;

	double recow = (xreco_max-xreco_min)/70.;
	eff_bkg = effbkg; 
	eff_sig = effsig; 


	for(int j=0;j<100;j++){
		double xgenval = xgen_min+  j*(xgen_max-xgen_min)/100.;
		xgencopy->setVal(xgenval);
		integral_store[j]=0.;
		bkgcont[j]= dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy);
		coscont[j]= cosspline.Eval(xgenval);
		sincont[j]= sinspline.Eval(xgenval);
		eff_bkgcont[j]= eff_bkg.Eval(xgenval);
		eff_sigcont[j]= eff_sig.Eval(xgenval);
		for(int i=0;i<71;i++){
			xrecocopy->setVal(xreco_min+  i*(xreco_max-xreco_min)/70.);
			Double_t gen = dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(*xrecocopy);
			resoval [i][j] = gen;
			integral_store[j] += gen*recow;
		}
		bins = xrecocopy->getBins();
	}


} 


Width_conv::Width_conv(const Width_conv& other, const char* name) :  
	RooAbsPdf(other,name), 
	xreco("xreco",this,other.xreco),
	mean("mean",this,other.mean),
	width("width",this,other.width),
	coupl("coupl",this,other.coupl),
	_coefList("coefList",this,other._coefList)

{ 
	xgen_min = other.xgen_min;
	xgen_max = other.xgen_max;
	xreco_min = other.xreco_min;
	xreco_max = other.xreco_max;
	for(int j=0;j<100;j++){
		bkgcont[j]=other.bkgcont[j];
		coscont[j]=other.coscont[j];
		sincont[j]=other.sincont[j];
		eff_bkgcont[j]=other.eff_bkgcont[j];
		eff_sigcont[j]=other.eff_sigcont[j];
		for(int i=0;i<71;i++){
			resoval[i][j]=other.resoval[i][j];
		}
		integral_store [j]= other.integral_store[j];
	}
	bins = other.bins;
	nbkg = other.nbkg;
	nsig = other.nsig;
	eff_bkg = other.eff_bkg; 
	eff_sig = other.eff_sig; 
	cosspline= other.cosspline;
	sinspline= other.sinspline;
} 


Double_t Width_conv::evaluate() const 
{ 
	double value = 0.;
	if(xreco>xreco_max)
		return 0;
	if(xreco<xreco_min)
		return 0;

	RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(mean.arg());
	TIterator *arglist= xlist->createIterator();
	RooRealVar *xgencopy= (RooRealVar*)arglist->Next();

	double genwid = (xgen_max-xgen_min)/100.;



	double left_edge = 0.;
	double left_width =  0;
	double right_edge = 0.;
	double right_width =  0;
//	int  re = (xreco-xreco_min)/((xreco_max-xreco_min)/70.); 

	{
		int nbinsC=40;
		double bwC= 0.5;
		if(width>1.){
			nbinsC=int(20*width*2);
			bwC=20./double(nbinsC);
		}
		double low_edge = mean-10*width-0.5*bwC*width;
		double high_edge = mean+10*width+0.5*bwC*width;
		int nl = (low_edge- xgen_min-0.5*genwid)/genwid;
		int nr = (high_edge-xgen_min-0.5*genwid)/genwid+1; 
		left_edge = nl*genwid+xgen_min+0.5*genwid;
		right_edge = nr*genwid+xgen_min+0.5*genwid;
		double left_center = (left_edge + low_edge )/2.;
		double right_center = (right_edge + high_edge )/2.;
		left_width = low_edge - left_edge;
		right_width = right_edge - high_edge;


		xgencopy->setVal(mean);
		//Double_t reso_s = dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(xreco.arg());


		//		double eff_T2 = eff_bkg.Eval(mean);

		for(int k=0;k<nbinsC+1;k++){
			double cvalue = mean-10*width + k*width*bwC;
			if(cvalue < xgen_min || cvalue>=xgen_max)
				continue;
			int  j = (cvalue-xgen_min)/genwid; 
			Double_t T4 = sincont[j];
			Double_t T3 = coscont[j];
			Double_t T2 = bkgcont[j];
			double eff_T2 = eff_bkgcont[j];
			double eff_T1 = eff_sigcont[j];
			double a= (cvalue*cvalue-mean*mean);
			double b = mean*width;
			double cossig = a/TMath::Sqrt(a*a+b*b);
			double sinsig = b/TMath::Sqrt(a*a+b*b);
			xgencopy->setVal( cvalue);
			Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
//			Double_t T2 = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy);
			double gen_s= nbkg*T2* eff_T2+ fabs(coupl) * nsig*T1*eff_T1 + 2.*sqrt(T2*T1*nbkg*nsig*coupl*eff_T2*eff_T1)*(cossig*T3+T4*sinsig); 
			Double_t reso_s = dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(xreco.arg());
			value+= reso_s*gen_s*bwC*width;
		}

		if(left_center> xgen_min && left_center <xgen_max){
		xgencopy->setVal(left_center);
		int  jleft = (left_center-xgen_min)/genwid; 
		Double_t reso_left= dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(xreco.arg());

		double eff_T2_left = eff_bkgcont[jleft];
		double eff_T1_left= eff_sigcont[jleft]; 
		Double_t T1_left = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
		Double_t T2_left = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy);
		Double_t T4_left = sincont[jleft];
		Double_t T3_left = coscont[jleft];
		double a_left= (left_center*left_center-mean*mean);
		double b_left = mean*width;
		double cossig_left = a_left/TMath::Sqrt(a_left*a_left+b_left*b_left);
		double sinsig_left = b_left/TMath::Sqrt(a_left*a_left+b_left*b_left);

		double gen_s_left= nbkg*T2_left*eff_T2_left+ fabs(coupl) * nsig*T1_left * eff_T1_left+ 2.*sqrt(eff_T2_left*T2_left*T1_left*eff_T1_left*nbkg*nsig*coupl)*(T3_left*cossig_left+T4_left*sinsig_left); 
		value+= reso_left*gen_s_left*left_width;
		}

		if(right_center> xgen_min && right_center <xgen_max){
		xgencopy->setVal(right_center);
		Double_t reso_right= dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(xreco.arg());
		int  jright = (right_center-xgen_min)/genwid; 
		double eff_T2_right = eff_bkgcont[jright];
		double eff_T1_right = eff_sigcont[jright]; 
		Double_t T1_right = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
		Double_t T2_right = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy); 
		Double_t T4_right = sincont[jright];
		Double_t T3_right = coscont[jright];
		double a_right= (right_center*right_center-mean*mean);
		double b_right = mean*width;
		double cossig_right = a_right/TMath::Sqrt(a_right*a_right+b_right*b_right);
		double sinsig_right = b_right/TMath::Sqrt(a_right*a_right+b_right*b_right);

		double gen_s_right= nbkg*T2_right*eff_T2_right+ fabs(coupl) * nsig*T1_right*eff_T1_right + 2.*sqrt(eff_T1_right*T2_right*T1_right*nbkg*eff_T2_right*nsig*coupl)*(T3_right*cossig_right+T4_right*sinsig_right); 
		//		double gen_s_right= nbkg*T2_right+ fabs(coupl) * nsig*T1_right + 2.*sqrt(T2_right*T1_right*nbkg*nsig*coupl)*(T3_right*cossig_right+T4_right*sinsig_right); 
		value+= reso_right*gen_s_right*right_width;
		}

	}

	for(int i=0;i<100;i++){
		double genval = xgen_min+  i*genwid;
		if(genval < xgen_min || genval>=xgen_max)
			continue;
		if(genval <= right_edge && genval >=left_edge){
			continue;
		}
		xgencopy->setVal(genval);

		Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);

		Double_t T2 = bkgcont[i]; 
		Double_t T3 = coscont[i]; 
		Double_t T4 = sincont[i];
		double eff_T2= eff_bkgcont[i]; 
		double eff_T1=eff_sigcont[i];

		double a= (genval*genval-mean*mean);
		double b = mean*width;
		double cossig = a/TMath::Sqrt(a*a+b*b);
		double sinsig = b/TMath::Sqrt(a*a+b*b);

		double gen= nbkg*T2 * eff_T2+ fabs(coupl) * nsig*T1 *eff_T1+ 2.*sqrt(eff_T2*eff_T1*T2*T1*nbkg*nsig*coupl)*(cossig*T3+sinsig*T4); 
		//double reso = resoval[re][i]; 
		double reso = dynamic_cast<const RooAbsPdf*>(_coefList.at(2))->getVal(xreco.arg());
		value += reso*gen*genwid;
	}
	return value;

} 

Int_t Width_conv::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

	if (matchArgs(allVars,analVars,xreco)) return 4 ;
	//	else return 3;
	//  if (matchArgs(allVars,analVars,xgen)) return 3 ;

	return 0 ;

}

Double_t Width_conv::analyticalIntegral(Int_t code, const char* rangeName) const
{
	switch(code)
	{

		case 4: 
			{
				double integral=0;
				RooArgSet *xlist= dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getParameters(mean.arg());
				TIterator *arglist= xlist->createIterator();
				RooRealVar *xgencopy= (RooRealVar*)arglist->Next();

				double genwid = (xgen_max-xgen_min)/100.;
		//		double low_edge = mean-10*width;
		//		double high_edge = mean+10*width;
				double left_edge = 0.;
				double left_width =  0;
				double right_edge = 0.;
				double right_width =  0;

				{
		int nbinsC=40;
		double bwC= 0.5;
	//	if(width>0.2){
	//		nbinsC=80;
	//		bwC=0.25;
	//	}
		if(width>1.){
			nbinsC=int(20*width*2);
			bwC=20./double(nbinsC);
		}
		double low_edge = mean-10*width-0.5*bwC*width;
		double high_edge = mean+10*width+0.5*bwC*width;
					int nl = (low_edge- xgen_min-0.5*genwid)/genwid;
					int nr = (high_edge-xgen_min-0.5*genwid)/genwid+1; 
					left_edge = nl*genwid+xgen_min+0.5*genwid;
					right_edge = nr*genwid+xgen_min+0.5*genwid;
					double left_center = (left_edge + low_edge )/2.;
					double right_center = (right_edge + high_edge )/2.;
					left_width = low_edge - left_edge;
					right_width = right_edge - high_edge;


		for(int k=0;k<nbinsC+1;k++){
						double cvalue = mean-10*width + k*width*bwC;
						if(cvalue < xgen_min || cvalue>=xgen_max)
							continue;
						int  j = (cvalue-xgen_min)/genwid; 
						Double_t T2 = bkgcont[j]; 
						Double_t T4 = sincont[j];
						Double_t T3 = coscont[j];
						double eff_T2 = eff_bkgcont[j];
						double eff_T1 = eff_sigcont[j];
						double a= (cvalue*cvalue-mean*mean);
						double b = mean*width;
						double cossig = a/TMath::Sqrt(a*a+b*b);
						double sinsig = b/TMath::Sqrt(a*a+b*b);
						xgencopy->setVal( cvalue);
						Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
						//Double_t T2 = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy);
						double gen_s= eff_T2*nbkg*T2 + fabs(coupl) * nsig*T1 * eff_T1+ 2.*sqrt(eff_T2*eff_T1*T2*T1*nbkg*nsig*coupl)*(cossig*T3+sinsig*T4); 
						integral+= gen_s*bwC*width*integral_store[j];
					}

					if(left_center> xgen_min && left_center <xgen_max){
					xgencopy->setVal(left_center);
					int  jleft = (left_center-xgen_min)/genwid; 
					double eff_T2_left = eff_bkgcont[jleft];
					double eff_T1_left= eff_sigcont[jleft]; 
					Double_t T2_left = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy); 
					Double_t T1_left = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
					Double_t T4_left = sincont[jleft];
					Double_t T3_left = coscont[jleft];
					double a_left= (left_center*left_center-mean*mean);
					double b_left = mean*width;
					double cossig_left = a_left/TMath::Sqrt(a_left*a_left+b_left*b_left);
					double sinsig_left = b_left/TMath::Sqrt(a_left*a_left+b_left*b_left);

					double gen_s_left= nbkg*T2_left*eff_T2_left+ fabs(coupl) * nsig*T1_left * eff_T1_left+ 2.*sqrt(eff_T2_left*T2_left*T1_left*eff_T1_left*nbkg*nsig*coupl)*(T3_left*cossig_left+T4_left*sinsig_left); 
					integral += gen_s_left*left_width*integral_store[jleft];
					}
					if(right_center> xgen_min && right_center <xgen_max){
					xgencopy->setVal(right_center);
					int  jright = (right_center-xgen_min)/genwid; 
					double eff_T2_right = eff_bkgcont[jright];
					double eff_T1_right = eff_sigcont[jright]; 
					Double_t T1_right = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
					Double_t T2_right = dynamic_cast<const RooAbsPdf*>(_coefList.at(1))->getVal(*xgencopy);
					Double_t T4_right = sincont[jright];
					Double_t T3_right = coscont[jright];
					double a_right= (right_center*right_center-mean*mean);
					double b_right = mean*width;
					double cossig_right = a_right/TMath::Sqrt(a_right*a_right+b_right*b_right);
					double sinsig_right = b_right/TMath::Sqrt(a_right*a_right+b_right*b_right);

					double gen_s_right= nbkg*T2_right*eff_T2_right+ fabs(coupl) * nsig*T1_right*eff_T1_right + 2.*sqrt(eff_T1_right*T2_right*T1_right*nbkg*eff_T2_right*nsig*coupl)*(T3_right*cossig_right+T4_right*sinsig_right); 
					integral+= gen_s_right*right_width*integral_store[jright];
					}

				}
				for(int i=0;i<100;i++){
					double genval = xgen_min+  i*genwid;
					if(genval < xgen_min || genval>=xgen_max)
						continue;
					if(genval <=right_edge && genval >=left_edge){
						continue;
					}
					xgencopy->setVal(genval);

					Double_t T1 = dynamic_cast<const RooAbsPdf*>(_coefList.at(0))->getVal(*xgencopy);
					Double_t T2 = bkgcont[i]; 
					Double_t T3 = coscont[i]; 
					Double_t T4 = sincont[i];
					double eff_T2= eff_bkgcont[i]; 
					double eff_T1=eff_sigcont[i];

					double a= (genval*genval-mean*mean);
					double b = mean*width;
					double cossig = a/TMath::Sqrt(a*a+b*b);
					double sinsig = b/TMath::Sqrt(a*a+b*b);

					double gen= nbkg*T2 * eff_T2+ fabs(coupl) * nsig*T1 *eff_T1+ 2.*sqrt(eff_T2*eff_T1*T2*T1*nbkg*nsig*coupl)*(cossig*T3+sinsig*T4); 
					integral += gen*genwid*integral_store[i];
				}


				return integral;
			}

	}

	assert(0) ;
	return 0 ;
}
