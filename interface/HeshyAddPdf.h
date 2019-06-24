#ifndef HiggsAnalysis_CombinedLimit_HeshyAddPdf_h
#define HiggsAnalysis_CombinedLimit_HeshyAddPdf_h

#include "RooAddPdf.h"

class HeshyAddPdf : public RooAddPdf {
public:
  HeshyAddPdf();
  HeshyAddPdf(const char *name, const char *title=0);

  HeshyAddPdf(const char *name, const char *title,
	    RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooAbsReal& coef1);
  HeshyAddPdf(const char *name, const char *title, const RooArgList& pdfList);
  HeshyAddPdf(const char *name, const char *title, const RooArgList& pdfList, const RooArgList& coefList, Bool_t recursiveFraction=kFALSE);
  
  HeshyAddPdf(const HeshyAddPdf& other, const char* name=0);
  virtual TObject* clone(const char* newname) const {return new HeshyAddPdf(*this,newname);}

private:
  ClassDef(HeshyAddPdf, 1) // PDF representing a sum of PDFs
};

#endif
