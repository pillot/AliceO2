// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHPreClustering/PreClusterFinderMapping.h"

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

std::vector<std::unique_ptr<Mapping::MpDE>> Mapping::readMapping(const char* mapFile)
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

    inFile.read(reinterpret_cast<char*>(&de.iCath[0]), sizeof(de.iCath[0]) * 2);
    inFile.read(reinterpret_cast<char*>(&de.nPads[0]), sizeof(de.nPads[0]) * 2);

    int nPadsInDE = de.nPads[0] + de.nPads[1];

    de.pads = std::make_unique<MpPad[]>(nPadsInDE);

    for (int ip = 0; ip < nPadsInDE; ++ip) {
      inFile.read(reinterpret_cast<char*>(&(de.pads[ip])), sizeof(de.pads[ip]));
      ++totalNumberOfPads;
    }

    int mapsize(2 * nPadsInDE);
    auto themap = std::make_unique<int64_t[]>(mapsize);
    int nMapElements(0);

    inFile.read(reinterpret_cast<char*>(themap.get()), sizeof(int64_t) * mapsize);

    for (int iPlane = 0; iPlane < 2; ++iPlane) {
      for (int ip = 0; ip < de.nPads[iPlane]; ++ip) {

        // Decode the unique ID value
        int cathode = (themap[nMapElements] & 0x40000000) >> 30;
        int deid = themap[nMapElements] & 0xFFF;
        int index = (themap[nMapElements] & 0x3FFFFFFF) >> 12;
        int manuid = index & 0xFFF;
        int manuch = (index >> 12) & 0x3F;

        // clear the cathode ID field in the unique ID value
        themap[nMapElements] &= 0x3FFFFFFF;
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
  static int manu2ds[64]={62,61,63,60,59,55,58,57,56,54,50,46,42,39,37,41,
      35,36,33,34,32,38,43,40,45,44,47,48,49,52,51,53,
      7, 6, 5, 4, 2, 3, 1, 0, 9,11,13,15,17,19,21,23,
      31,30,29,28,27,26,25,24,22,20,18,16,14,12,10, 8};

  int NDES = 228;
  int DEvec[228] = {100, 101, 102, 103, 200, 201, 202, 203, 300, 301, 302, 303, 400, 401, 402, 403, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 609, 610, 611, 612, 613, 614, 615, 616, 617, 600, 601, 602, 603, 604, 605, 606, 607, 608, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709, 710, 711, 712, 713, 714, 715, 716, 717, 718, 719, 720, 721, 722, 723, 724, 725, 800, 801, 802, 803, 804, 805, 806, 807, 808, 809, 810, 811, 812, 813, 814, 815, 816, 817, 818, 819, 820, 821, 822, 823, 824, 825, 907, 908, 909, 910, 911, 912, 913, 914, 915, 916, 917, 918, 919, 920, 921, 922, 923, 924, 925, 900, 901, 902, 903, 904, 905, 906, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308, 1309, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1400, 1401, 1402, 1403, 1404, 1405, 1406, 1407, 1408, 1409, 1410, 1411, 1412, 1413, 1414, 1415, 1416, 1417};

  std::vector<std::unique_ptr<Mapping::MpDE>> detectionElements{};

  int totalNumberOfPads(0);

  int numberOfDetectionElements(0);

  // Loop over all possible detection element IDs
  // The non-existing ones rise an exception and are be skipped
  for(int deidx = 0; deidx < NDES; deidx ++) {

    int deid = DEvec[deidx];
    try {
      // Get the Segmentation object for the current detection element IDs
      // The constructor rises an exception if the deid is unknown
      mapping::Segmentation segment(deid);
      const mapping::CathodeSegmentation& bend = segment.bending();
      const mapping::CathodeSegmentation& nbend = segment.nonBending();

      detectionElements.push_back(std::make_unique<Mapping::MpDE>());
      Mapping::MpDE& de( *(detectionElements.back()) );

      de.uid = deid;
      de.segment = std::make_unique<mapping::Segmentation>(deid);
      de.iCath[0] = 0;
      de.iCath[1] = 1;

      de.nPads[0] = bend.nofPads();
      de.nPads[1] = nbend.nofPads();

      int nPadsInDE = de.nPads[0] + de.nPads[1];

      de.pads = std::make_unique<MpPad[]>(nPadsInDE);
      std::unique_ptr<uint32_t[]> pad_ids = std::make_unique<uint32_t[]>(nPadsInDE);

      // loop on the pads of the current segmentation, and insert each of them into the mapping
      int ip{0};
      auto addPadId = [&readoutVersion,&segment,&de,&pad_ids,&ip,&totalNumberOfPads](int padid) {

        // set pad coordinates
        double padX = segment.padPositionX(padid);
        double padY = segment.padPositionY(padid);
        double padSizeX = segment.padSizeX(padid);
        double padSizeY = segment.padSizeY(padid);
        de.pads[ip].area[0][0] = padX - padSizeX/2;
        de.pads[ip].area[0][1] = padX + padSizeX/2;
        de.pads[ip].area[1][0] = padY - padSizeY/2;
        de.pads[ip].area[1][1] = padY + padSizeY/2;

        pad_ids[ip] = padid;

        // build the unique ID of the pad
        bool isBending = segment.isBendingPad(padid);
        int iCath = (isBending ? 0 : 1);
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
        uid += (manuid << 12) & 0xFFF000;
        uid += (manuch << 24) & 0x3F000000;
        //uid += (iCath << 30) & 0x40000000;

        // insert the pad into the hash table
        de.padIndices[iPlane].Add(uid, ip+1);

        ++ip;
        ++totalNumberOfPads;
      };
      segment.forEachPad(addPadId);

      // fill vector of pad neighbours
      for(ip = 0; ip < nPadsInDE; ip++) {

        de.pads[ip].nNeighbours = 0;
        uint32_t cpadid = pad_ids[ip];
        bool isBending = segment.isBendingPad(cpadid);
        int iPlane = (isBending ? 0 : 1);

        auto addPadNeighbours = [&segment,&cpadid,&iPlane,&de,&ip](int padid) {

          if( de.pads[ip].nNeighbours == 10 ) return;

          // build the unique ID of the pad
          bool isBending = segment.isBendingPad(padid);
          int iCath = (isBending ? 0 : 1);
          int iPlane = (isBending ? 0 : 1);
          int deid = segment.detElemId();
          int manuid = segment.padDualSampaId(padid);
          int manuch = segment.padDualSampaChannel(padid);
          uint32_t uid = deid & 0xFFF;
          uid += (manuid << 12) & 0xFFF000;
          uid += (manuch << 24) & 0x3F000000;
          //uid += (iCath << 30) & 0x40000000;

          auto iPad = de.padIndices[iPlane].GetValue(uid) - 1;
          de.pads[ip].neighbours[ de.pads[ip].nNeighbours ] = iPad;
          de.pads[ip].nNeighbours++;
        };
        segment.forEachNeighbouringPad(cpadid, addPadNeighbours);
      }

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
