// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
///  HAND CODED
///

#include "CH.cxx"
void fillElec2DetCH8L(std::map<uint16_t, uint32_t>& e2d)
{
  add(e2d, 819, 108, 861, 0, 0);
  add(e2d, 819, 107, 861, 0, 2);
  add(e2d, 819, 106, 861, 0, 4);
  add(e2d, 819, 1133, 861, 1, 0);
  add(e2d, 819, 1134, 861, 1, 2);
}

void fillSolar2CruLinkCH8L(std::map<uint16_t, uint32_t>& s2c)
{
  add_cru(s2c, 0, 0, 861, 819);
}

