// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHBase/Mapping.h"
#include "MCHMappingInterface/Segmentation.h"

#include <cassert>
#include <fstream>
#include <iostream>

#include <TMath.h>

#include <FairMQLogger.h>

namespace o2
{
namespace mch
{

using namespace std;

static int manu2ds[64]={62,61,63,60,59,55,58,57,56,54,50,46,42,39,37,41,
    35,36,33,34,32,38,43,40,45,44,47,48,49,52,51,53,
    7, 6, 5, 4, 2, 3, 1, 0, 9,11,13,15,17,19,21,23,
    31,30,29,28,27,26,25,24,22,20,18,16,14,12,10, 8};



std::vector<std::unique_ptr<Mapping::MpDE>> Mapping::readMapping(const char* mapFile, int readoutVersion)
{
  std::vector<std::unique_ptr<Mapping::MpDE>> detectionElements{};

  ifstream inFile(mapFile, ios::binary);

  if (!inFile.is_open()) {
    LOG(ERROR) << "Can not open " << mapFile;
    return detectionElements;
  }

  int totalNumberOfPads(0);

  int numberOfDetectionElements(0);
  inFile.read(reinterpret_cast<char*>(&numberOfDetectionElements), sizeof(numberOfDetectionElements));

  LOG(INFO) << "numberOfDetectionElements = " << numberOfDetectionElements;

  detectionElements.reserve(numberOfDetectionElements);

  for (int i = 0; i < numberOfDetectionElements; ++i) {

    detectionElements.push_back(std::make_unique<Mapping::MpDE>());
    Mapping::MpDE& de(*detectionElements[i]);

    inFile.read(reinterpret_cast<char*>(&de.uid), sizeof(de.uid));

    bool verb = (de.uid==819) ? true : false;

    inFile.read(reinterpret_cast<char*>(&de.iCath[0]), sizeof(de.iCath[0]) * 2);
    inFile.read(reinterpret_cast<char*>(&de.nPads[0]), sizeof(de.nPads[0]) * 2);

    if(verb) {
      printf("[readMapping] cathode indexes: %d,%d\n", (int)de.iCath[0], (int)de.iCath[1]);
      printf("[readMapping] number of pads: %d,%d\n", (int)de.nPads[0], (int)de.nPads[1]);
    }

    int nPadsInDE = de.nPads[0] + de.nPads[1];

    de.pads = std::make_unique<MpPad[]>(nPadsInDE);

    for (int ip = 0; ip < nPadsInDE; ++ip) {
      inFile.read(reinterpret_cast<char*>(&(de.pads[ip])), sizeof(de.pads[ip]));
      ++totalNumberOfPads;
      if(verb) {
        printf("[readMapping] de.pads[%d]: %f,%f -> %f,%f\n", ip,
            de.pads[ip].area[0][0],de.pads[ip].area[1][0],
            de.pads[ip].area[0][1],de.pads[ip].area[1][1]);
      }
    }

    int mapsize(2 * nPadsInDE);
    auto themap = std::make_unique<int64_t[]>(mapsize);
    int nMapElements(0);

    inFile.read(reinterpret_cast<char*>(themap.get()), sizeof(int64_t) * mapsize);

    for (int iPlane = 0; iPlane < 2; ++iPlane) {
      for (int ip = 0; ip < de.nPads[iPlane]; ++ip) {
        if( readoutVersion == 2 ) {

          if(verb) printf("[readMapping] de.padIndices[%d].Add(%lu,%lu) \n", iPlane,
              themap[nMapElements],themap[nMapElements+1]);

          // Decode the unique ID value
          int cathode = (themap[nMapElements] & 0x40000000) >> 30;
          int deid = themap[nMapElements] & 0xFFF;
          int index = (themap[nMapElements] & 0x3FFFFFFF) >> 12;
          int manuid = index & 0xFFF;
          int manuch = (index >> 12) & 0x3F;
          if(verb) printf("              DEid=%d cathode=%d index=%d manuid=%d manuch=%d\n",
              deid, cathode, index, manuid, manuch);

          // DualSAMPA readout boards need a re-mapping of the channel numbers
          int dsch = manu2ds[manuch];
          // reset the channel number portion of the unique ID
          themap[nMapElements] &= 0x40FFFFFF;
          // set the channel number according to the DualSAMPA mapping
          themap[nMapElements] |= ( (dsch << 24) & 0x3F000000 );
          if(verb) printf("              new unique ID after DS re-mapping (%d -> %d): %d\n",
              manuch, dsch, (int)themap[nMapElements]);

          if(verb) {
            int ip = themap[nMapElements + 1] - 1;
            printf("              de.pads[%d]: %f,%f -> %f,%f\n", ip,
                de.pads[ip].area[0][0],de.pads[ip].area[1][0],
                de.pads[ip].area[0][1],de.pads[ip].area[1][1]);
          }
        }
        de.padIndices[iPlane].Add(themap[nMapElements], themap[nMapElements + 1]);
        nMapElements += 2;
      }
    }

    assert(nMapElements == 2 * nPadsInDE);
    assert(de.padIndices[0].GetSize() + de.padIndices[1].GetSize() == nPadsInDE);
  }

  if (totalNumberOfPads != 1064008 + 20992) {
    LOG(ERROR) << "totalNumberOfPads = " << totalNumberOfPads << "!= from the expected " << 1064008 + 20992;
    detectionElements.clear();
  }

  return detectionElements;
}


//_________________________________________________________________________________________________
std::vector<std::unique_ptr<Mapping::MpDE>> Mapping::loadO2Mapping(int readoutVersion)
{
  std::vector<std::unique_ptr<Mapping::MpDE>> detectionElements{};

  int totalNumberOfPads(0);

  int numberOfDetectionElements(0);

  // Loop over all possible detection element IDs
  // The non-existing ones rise an exception and are be skipped
  for(int deid = 0; deid <= 5000; deid ++) {

    try {
      // Get the Segmentation object for the current detection element IDs
      // The constructor rises an exception if the deid is unknown
      mapping::Segmentation segment(deid);
      const mapping::CathodeSegmentation& bend = segment.bending();
      const mapping::CathodeSegmentation& nbend = segment.nonBending();

      detectionElements.push_back(std::make_unique<Mapping::MpDE>());
      Mapping::MpDE& de( *(detectionElements.back()) );

      printf("[readMappingO2] loading DE %d\n", deid);

      de.uid = deid;
      de.iCath[0] = 1;
      de.iCath[1] = 0;

      de.nPads[0] = bend.nofPads();
      de.nPads[1] = nbend.nofPads();

      int nPadsInDE = de.nPads[0] + de.nPads[1];

      de.pads = std::make_unique<MpPad[]>(nPadsInDE);

      // loop on the pads of the current segmentation, and insert each of them into the mapping
      int ip{0};
      auto addPadId = [&readoutVersion,&segment,&de,&ip,&totalNumberOfPads](int padid) {

        // set pad coordinates
        float padX = segment.padPositionX(padid);
        float padY = segment.padPositionY(padid);
        float padSizeX = segment.padSizeX(padid);
        float padSizeY = segment.padSizeY(padid);
        de.pads[ip].area[0][0] = padX - padSizeX/2;
        de.pads[ip].area[0][1] = padX + padSizeX/2;
        de.pads[ip].area[1][0] = padY - padSizeY/2;
        de.pads[ip].area[1][1] = padY + padSizeY/2;
        //printf("[readMappingO2] de.pads[%d]: %f,%f -> %f,%f\n", ip,
        //       de.pads[ip].area[0][0],de.pads[ip].area[1][0],
        //       de.pads[ip].area[0][1],de.pads[ip].area[1][1]);

        // build the unique ID of the pad
        bool isBending = segment.isBendingPad(padid);
        int iCath = (isBending ? 1 : 0);
        int iPlane = (isBending ? 0 : 1);
        int deid = segment.detElemId();
        int manuid = segment.padDualSampaId(padid);
        int manuch = segment.padDualSampaChannel(padid);
        if( readoutVersion == 2 ) {
          // DualSAMPA readout boards need a re-mapping of the channel numbers
          int dsch = manu2ds[manuch];
          manuch = dsch;
        }

        uint32_t uid = deid & 0xFFF;
        //uid += (padid << 12) & 0x3FFFF000;
        uid += (manuid << 12) & 0xFFF000;
        uid += (manuch << 24) & 0x3F000000;
        uid += (iCath << 30) & 0x40000000;

        // insert the pad into the hash table
        de.padIndices[iPlane].Add(uid, padid);
        //printf("                isBending=%d  iPlane=%d  iCath=%d  manuid=%d  manuch=%d  uid=%u\n",
        //       (int)isBending, iPlane, iCath, manuid, manuch, uid);

        ++ip;
        ++totalNumberOfPads;
      };
      segment.forEachPad(addPadId);

      // TODO: fill vector of pad neighbours

      //assert(de.padIndices[0].GetSize() + de.padIndices[1].GetSize() == nPadsInDE);

    }
    catch(std::exception& e) {
      continue;
    }
  }

  if (false && totalNumberOfPads != (1064008 + 20992)) {
    LOG(ERROR) << "totalNumberOfPads = " << totalNumberOfPads << "!= from the expected " << 1064008 + 20992;
    detectionElements.clear();
  }

  return detectionElements;
}

//_________________________________________________________________________________________________
bool Mapping::areOverlapping(float area1[2][2], float area2[2][2], float precision)
{
  /// check if the two areas overlap
  /// precision in cm: positive = increase pad size / negative = decrease pad size

  if (area1[0][0] - area2[0][1] > precision) {
    return false;
  }
  if (area2[0][0] - area1[0][1] > precision) {
    return false;
  }
  if (area1[1][0] - area2[1][1] > precision) {
    return false;
  }
  if (area2[1][0] - area1[1][1] > precision) {
    return false;
  }

  return true;
}

//_________________________________________________________________________________________________
bool Mapping::areOverlappingExcludeCorners(float area1[2][2], float area2[2][2])
{
  /// check if the two areas overlap (excluding pad corners)

  // precision in cm: positive = increase pad size / negative = decrease pad size
  constexpr float precision = 1.e-4;

  if (areOverlapping(area1, area2, precision)) {
    for (int ip1 = 0; ip1 < 2; ++ip1) {
      for (int ip2 = 0; ip2 < 2; ++ip2) {
        if (TMath::Abs(area1[0][ip1] - area2[0][1 - ip1]) < precision &&
            TMath::Abs(area1[1][ip2] - area2[1][1 - ip2]) < precision) {
          return false;
        }
      }
    }
    return true;
  }

  return false;
}

} // namespace mch
} // namespace o2
