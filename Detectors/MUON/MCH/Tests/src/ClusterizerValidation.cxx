// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <fstream>
#include <cstdio>

#include "MCHBase/Digit.h"
//#include "MCHBase/DigitBlock.h"
#include "MCHPreClustering/PreClusterBlock.h"
#include "MCHPreClustering/PreClusterFinder.h"
#include "DigitsFileReader.h"
#include "MCHClustering/ClusteringForTest.h"
#include "Validation.h"

#include "TMath.h"
#include "TRandom.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include <cmath>
#include <TApplication.h>


using namespace o2::mch;
using namespace std;

int main(int argc, char** argv){
    
    Validation validation;
    
    TApplication app ("app",&argc,argv);
    cout << "\n\n==========\nRunning the Validation procedure of (pre)clustering" << endl;
    
    cout << "\n\n==========\nGetting info for Bending plane\n\n" << endl;
   validation.InfoDE819b();
    cout << "\n\n==========\nGetting info for Non-Bending plane\n\n" << endl;
   validation.InfoDE819nb();
    cout << "\n\n==========\nHit generation, histograms plotting and digitization\n\n" << endl;
   validation.PlotMathieson2D();
    cout << "\n\n==========\nTesting the (pre)clustering\n\n" << endl;
   validation.TestClustering();
    
    cout << "\n\n==========\nValidation procedure terminated\n\n" << endl;
    
   //ResidualsCOG();
    
   // ResidualsCompare();
    
   app.Run(kTRUE);
    
    return 0;
}

