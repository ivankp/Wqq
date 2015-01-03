#ifndef TMinuitObj_h
#define TMinuitObj_h

#include "TMinuit.h"

class MnFcnBase {
  protected:
    Int_t npar;
    Double_t *par, *err;

  public:
    MnFcnBase();
    virtual ~MnFcnBase();

    Double_t get_par(Int_t i) const;
    Double_t get_err(Int_t i) const;

    virtual void operator()
    (Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) const =0;

  friend class TMinuitObj;
};

class TMinuitObj : public TMinuit {
  private:
    MnFcnBase* fFCNobj;

  public:
    TMinuitObj(Int_t maxpar);
    virtual ~TMinuitObj();

    virtual void SetFCN(MnFcnBase* fcn);
    virtual void SetFCN(void *fcn) { TMinuit::SetFCN(fcn); }
    virtual void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &f, Double_t *, Int_t))
    { TMinuit::SetFCN(fcn); }

    Int_t Eval(Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

    virtual Int_t Migrad();

};

#endif
