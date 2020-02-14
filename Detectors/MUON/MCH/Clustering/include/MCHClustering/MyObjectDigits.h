#ifndef ALICEO2_MCH_MYOBJECTDIGITS_H_
#define ALICEO2_MCH_MYOBJECTDIGITS_H_

#include "TObject.h"
#include <iostream>
#include <vector>
#include "MCHBase/Digit.h"


namespace o2
{
namespace mch
{
    
using namespace std;

//class MyObjectDigits : public TObject
//{
//
//public:
//
//    MyObjectDigits() = default;
//    virtual ~MyObjectDigits() = default;
//
//    std::vector<Digit> fDigits;
//    void setfDigits(std::vector<Digit> &digits) { fDigits = digits; }
//
//
//    ClassDef(MyObjectDigits,1); // A value of 1 indicates that this class is streamable.
//
//};

class MyObjectDigits : public TObject
{
 public:
  MyObjectDigits() {};
  ~MyObjectDigits() {};

  std::vector<Digit> getfDigits() const { return fDigits; }
  void setfDigits(std::vector<Digit> vectdigits) { fDigits = vectdigits; }

 private:
    std::vector<Digit> fDigits;

    ClassDef(MyObjectDigits, 1);
};

}
}

#endif // ALICEO2_MCH_MYOBJECTDIGITS_H_
