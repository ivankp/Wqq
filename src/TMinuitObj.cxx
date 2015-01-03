#include "TMinuitObj.h"

TMinuitObj::TMinuitObj(Int_t maxpar) // constructor
: TMinuit(maxpar) { }

TMinuitObj::~TMinuitObj() // destructor
{ }

void TMinuitObj::SetFCN(MnFcnBase* fcn) { fFCNobj = fcn; }

Int_t TMinuitObj::Eval(Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
  if (fFCNobj) fFCNobj->operator()(npar,grad,fval,par,flag);
  else printf("!!! Error in TMinuitObj::Eval() : fFCNobj pointer is null! \n");
  return 0;
}

MnFcnBase::MnFcnBase() : par(NULL), err(NULL) { }

MnFcnBase::~MnFcnBase()
{
  if (par) delete[] par;
  if (err) delete[] err;
}

Int_t TMinuitObj::Migrad() {
  const Int_t status = TMinuit::Migrad();
  fFCNobj->npar = GetNumPars();
  fFCNobj->par = new Double_t[fFCNobj->npar];
  fFCNobj->err = new Double_t[fFCNobj->npar];
  for (Int_t i=0;i<fFCNobj->npar;i++)
    GetParameter(i,fFCNobj->par[i],fFCNobj->err[i]);
  return status;
}


Double_t MnFcnBase::get_par(Int_t i) const {
  if (i<npar) return par[i];
  else return 0.;
}
Double_t MnFcnBase::get_err(Int_t i) const {
  if (i<npar) return err[i];
  else return 0.;
}
