#include "HiggsAnalysis/CombinedLimit/interface/HeshyAddPdf.h"

HeshyAddPdf::HeshyAddPdf() :
  RooAddPdf()
{
  delete[] _coefCache;
  _coefCache = new Double_t[1000];
}

HeshyAddPdf::HeshyAddPdf(const char *name, const char *title) :
  RooAddPdf(name, title)
{
  delete[] _coefCache;
  _coefCache = new Double_t[1000];
}

HeshyAddPdf::HeshyAddPdf(const char *name, const char *title,
            RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooAbsReal& coef1) :
  RooAddPdf(name, title, pdf1, pdf2, coef1)
{}

HeshyAddPdf::HeshyAddPdf(const char *name, const char *title, const RooArgList& pdfList) :
  RooAddPdf(name, title, pdfList)
{}

HeshyAddPdf::HeshyAddPdf(const char *name, const char *title, const RooArgList& pdfList, const RooArgList& coefList, Bool_t recursiveFraction) :
  RooAddPdf(name, title, pdfList, coefList, recursiveFraction)
{}

HeshyAddPdf::HeshyAddPdf(const HeshyAddPdf& other, const char* name) :
  RooAddPdf(other, name)
{}
