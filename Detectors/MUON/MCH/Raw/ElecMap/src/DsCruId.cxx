// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MCHRawElecMap/CruLinkId.h"
#include "MCHRawElecMap/DsCruId.h"
#include "Assertions.h"
#include <fmt/format.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

namespace o2::mch::raw
{
DsCruId::DsCruId(uint16_t linkId, uint8_t elinkGroupId, uint8_t elinkIndex)
  : mLinkId{linkId}, mElinkGroupId{elinkGroupId}, mElinkIndexInGroup{elinkIndex}
{
  impl::assertIsInRange("elinkGroupId", mElinkGroupId, 0, 7);
  impl::assertIsInRange("elinkIndex", mElinkIndexInGroup, 0, 4);
}

uint16_t encode(const DsCruId& id)
{
  return (id.linkId() & 0x3FF) | ((id.elinkGroupId() & 0x7) << 10) |
         ((id.elinkIndexInGroup() & 0x7) << 13);
}

std::optional<DsCruId> decodeDsCruId(uint16_t code)
{
  uint16_t linkId = code & 0x3FF;

  uint8_t groupId = (code & 0x1C00) >> 10;

  uint8_t index = (code & 0xE000) >> 13;

  if (groupId > 7) {
    return std::nullopt;
  }
  if (index > 4) {
    return std::nullopt;
  }
  return DsCruId(linkId, groupId, index);
}

std::optional<DsCruId> decodeDsCruId(std::string rep)
{
  std::istringstream is(rep);
  std::string line;
  std::vector<std::string> tokens;
  while (getline(is, line, '-')) {
    tokens.emplace_back(line);
  }
  if (tokens.size() < 3) {
    // not a valid representation of a DsCruId
    return std::nullopt;
  }
  if (tokens[0].empty() || tokens[0][0] != 'S') {
    // token[0] is not a valid representation of a linkId
    return std::nullopt;
  }
  if (tokens[1].empty() || tokens[1][0] != 'J') {
    // token[1] is not a valid representation of a groupId
    return std::nullopt;
  }
  if (tokens[2].size() < 3 || tokens[2][0] != 'D' || tokens[2][1] != 'S') {
    // token is not a valid representation of a DS
    return std::nullopt;
  }
  uint16_t linkId = std::atoi(tokens[0].substr(1).c_str());
  uint8_t groupId = std::atoi(tokens[1].substr(1).c_str());
  uint8_t index = std::atoi(tokens[2].substr(2).c_str());
  return DsCruId(linkId, groupId, index);
}

std::ostream& operator<<(std::ostream& os, const DsCruId& id)
{
  std::cout << fmt::format("DsCruId(SOLAR=S{:4d} GROUP=J{:2d} INDEX=DS{:2d}) CODE={:8d}",
                           id.linkId(), id.elinkGroupId(), id.elinkIndexInGroup(), encode(id));
  return os;
}

std::string asString(DsCruId dsId)
{
  CruLinkId cruLinkId(decodeCruLinkId(dsId.linkId()));
  return fmt::format("C{}-L{}-J{}-DS{}", cruLinkId.cruId(), cruLinkId.linkId(), dsId.elinkGroupId(), dsId.elinkIndexInGroup());
}

std::optional<uint8_t> groupFromElinkId(uint8_t elinkId)
{
  if (elinkId < 40) {
    return elinkId / 5;
  }
  return std::nullopt;
}

std::optional<uint8_t> indexFromElinkId(uint8_t elinkId)
{
  if (elinkId < 40) {
    return elinkId - (elinkId / 5) * 5;
  }
  return std::nullopt;
}
} // namespace o2::mch::raw
