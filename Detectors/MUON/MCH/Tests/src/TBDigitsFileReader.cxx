// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <iostream>

#include "TBDigitsFileReader.h"
#include "MCHBase/DigitBlock.h"
#include "MCHMappingInterface/Segmentation.h"

using namespace o2::mch;
using namespace std;


TBDigitsFileReader::TBDigitsFileReader(std::string inputFileName)
{
  mInputFile.open(inputFileName, ios::binary);
  if (!mInputFile.is_open()) {
    throw invalid_argument("Cannot open input file " + inputFileName);
  }
}


bool TBDigitsFileReader::readDigitsFromFile()
{
  int manu2ds[64]={62,61,63,60,59,55,58,57,56,54,50,46,42,39,37,41,
      35,36,33,34,32,38,43,40,45,44,47,48,49,52,51,53,
      7, 6, 5, 4, 2, 3, 1, 0, 9,11,13,15,17,19,21,23,
      31,30,29,28,27,26,25,24,22,20,18,16,14,12,10, 8};

  int ds2manu[64];
  for(int i = 0; i < 64; i++) {
    for(int j = 0; j < 64; j++) {
      if( manu2ds[j] != i ) continue;
      ds2manu[i] = j;
      break;
    }
  }


  digits.clear();

  /// send the digits of the current event

  int nDE = 0;
  int DE1 = 0; 
  mInputFile.read(reinterpret_cast<char*>(&nDE), sizeof(int));
  if(mInputFile.eof()) return false;

  for(int iDE = 0; iDE < nDE; iDE++) {

    int DE;
    mInputFile.read(reinterpret_cast<char*>(&DE), sizeof(int));
    float xtrk, ytrk;
    mInputFile.read(reinterpret_cast<char*>(&xtrk), sizeof(float));
    mInputFile.read(reinterpret_cast<char*>(&ytrk), sizeof(float));
    trkx[DE] = xtrk;
    trky[DE] = ytrk;

    int npad;
    mInputFile.read(reinterpret_cast<char*>(&npad), sizeof(int));

    int dsId, dsCh, size, time;
    float charge;
    for(int ih = 0; ih < npad; ih++) {

      mInputFile.read(reinterpret_cast<char*>(&dsId), sizeof(int));
      mInputFile.read(reinterpret_cast<char*>(&dsCh), sizeof(int));
      mInputFile.read(reinterpret_cast<char*>(&size), sizeof(int));
      mInputFile.read(reinterpret_cast<char*>(&time), sizeof(int));
      mInputFile.read(reinterpret_cast<char*>(&charge), sizeof(float));

      printf("B  hit %d  dsid=%d chan=%d  charge=%f  time=%d\n",
          ih, dsId, dsCh, charge, time);

      uint16_t adc = static_cast<uint16_t>(charge);

      uint16_t detId = DE;
      uint16_t dualSampaId = dsId;
      uint16_t dualSampaChannel = ds2manu[dsCh];

      try {
        mapping::Segmentation segment(detId);
        int padId = segment.findPadByFEE(dualSampaId, dualSampaChannel);
        if(padId < 0) continue;
        //digits.push_back( std::make_unique<Digit>(time, detId, padId, adc) );
        digits.push_back( std::make_unique<Digit>() );
        Digit* mchdigit = digits.back().get();
        mchdigit->setDetID(detId);
        mchdigit->setPadID(padId);
        mchdigit->setADC(adc);
        mchdigit->setTimeStamp(time);
      }
      catch(std::exception& e) {
        continue;
      }
    }
  }

  return true;
}


ssize_t TBDigitsFileReader::getNumberOfDigits()
{
  return digits.size();
}


void TBDigitsFileReader::storeDigits(void* bufferPtr)
{
  Digit* ptr = (Digit*)bufferPtr;
  for(unsigned int di = 0; di < digits.size(); di++) {

    memcpy(ptr, digits[di].get(), sizeof(Digit));
    ptr += 1;
  }
}
