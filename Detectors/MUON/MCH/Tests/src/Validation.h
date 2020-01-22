// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef Validation_h
#define Validation_h

#include <memory>
#include <fstream>
#include <stdio.h>

#include "MCHBase/Digit.h"

using namespace o2::mch;
using namespace std;

namespace o2 {

namespace mch {

 Double_t myMathieson2D(Double_t *x, Double_t *par);
Double_t myMathieson2D2hits(Double_t *x, Double_t *par);
 void myMath();
void ResidualsCOG();
void ResidualsCompare();

class Validation
{
public:
  Validation();
  void PlotMathieson2D();
  void InfoDE819b();
  void InfoDE819nb();
  void TestClustering();
    ssize_t getNumberOfDigits();
    void storeDigits(void* bufferPtr);
    
private:
    
    vector<double> lowxsb;
    vector<double> lowysb;
    vector<double> lowxsnb;
    vector<double> lowysnb;
    
    std::vector< std::unique_ptr<mch::Digit> > digits;
    mch::Digit* digitsBuffer;
    int nDigits;

};

}
}

#endif /* Validation_h */
