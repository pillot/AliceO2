#ifndef ALICEO2_MCH_MYOBJECTPARAMETERS_H_
#define ALICEO2_MCH_MYOBJECTPARAMETERS_H_

#include <TROOT.h>
#include <iostream>
#include <vector>
#include <TObject.h>

namespace o2
{
namespace mch
{

using namespace std;

//class MyObjectParameters : public TObject
//{
//
//public:
//
//    MyObjectParameters() = default;
//    MyObjectParameters(int cathode, double chargetot, double Ky3, double Kx3, int Nbdigits);
//    virtual ~MyObjectParameters() = default;
//
//    int fCathode;
//    Double_t fChargetot;
//    Double_t fKy3;
//    Double_t fKx3;
//    Int_t fNbdigits;
//
//    void setfCathode(int cathode) { fCathode = cathode; }
//    void setfChargetot(double chargetot) { fChargetot = chargetot; }
//    void setfKy3(double Ky3) { fKy3 = Ky3; }
//    void setfKx3(double Kx3) { fKx3 = Kx3; }
//    void setfNbdigits(int Nbdigits) { fNbdigits = Nbdigits; }
//
//
//
//
//    ClassDef(MyObjectParameters,1); // A value of 1 indicates that this class is streamable.
//
//
//};

class MyObjectParameters : public TObject
{
 public:
  MyObjectParameters() = default;
  MyObjectParameters(int cathode, double chargetot, double Ky3, double Kx3, int Nbdigits);
  virtual ~MyObjectParameters() = default;

  int getfCathode() const { return fCathode; }
  void setfCathode(int cathode) { fCathode = cathode; }
    double getfChargetot() const { return fChargetot; }
    void setfChargetot(double chargetot) { fChargetot = chargetot; }
    double getfKy3() const { return fKy3; }
    void setfKy3(double Ky3) { fKy3 = Ky3; }
    double getfKx3() const { return fKx3; }
    void setfKx3(double Kx3) { fKx3 = Kx3; }
    int getfNbdigits() const { return fNbdigits; }
    void setfNbdigits(int Nbdigits) { fNbdigits = Nbdigits; }

 private:
    int fCathode;
    double fChargetot;
    double fKy3;
    double fKx3;
    int fNbdigits;

  ClassDef(MyObjectParameters, 1);
};

}
}

#endif // ALICEO2_MCH_MYOBJECTPARAMETERS_H_
