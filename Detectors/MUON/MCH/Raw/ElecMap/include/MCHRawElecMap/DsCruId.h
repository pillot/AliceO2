// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_RAW_ELECMAP_DS_CRU_ID_H
#define O2_MCH_RAW_ELECMAP_DS_CRU_ID_H

#include <cstdint>
#include <iosfwd>
#include <optional>

namespace o2::mch::raw
{
/// A DsCruId is a tuple (cruLink,group,index) that identifies
/// a dual sampa from the CRU point of view.
class DsCruId
{
 public:
  explicit DsCruId(uint16_t linkId, uint8_t elinkGroupId, uint8_t elinkIndex);

  /// The eLinkIndexInGroup is a number between 0 and 4
  /// and represents the index of a dual sampa within a given flex (aka group)
  constexpr uint8_t elinkIndexInGroup() const
  {
    return mElinkIndexInGroup;
  }

  /// The elinkGroupId is a number between 0 and 7
  /// Note that on the detector cable J = elinkGroupId+1
  constexpr uint8_t elinkGroupId() const
  {
    return mElinkGroupId;
  }

  /// elinkId is between 0 and 39 and represents one of the 40 dual sampas
  /// connected to a Solar (aka GBT)
  constexpr uint8_t elinkId() const
  {
    return mElinkGroupId * 5 + mElinkIndexInGroup;
  }

  /// solarId is an identifier that uniquely identify a solar board
  constexpr uint16_t linkId() const
  {
    return mLinkId;
  }

  bool operator==(const DsCruId& rhs)
  {
    return mLinkId == rhs.mLinkId &&
           mElinkIndexInGroup == rhs.mElinkIndexInGroup &&
           mElinkGroupId == rhs.mElinkGroupId;
  }
  bool operator!=(const DsCruId& rhs)
  {
    return !(*this == rhs);
  }

 private:
  uint16_t mLinkId;           // [0..63] * [0..11]
  uint8_t mElinkGroupId;      // 0..7
  uint8_t mElinkIndexInGroup; // 0..4
};

uint16_t encode(const DsCruId& id);

/// Creates (if possible) a DsCruId object from a code.
/// If the code does not correspond to a valid DsElectId
/// std::nullopt is returned
std::optional<DsCruId> decodeDsCruId(uint16_t code);

/// Creates (if possible) a DsCruId object from a string representation
/// If the string does not correspond to a valid DsElectId
/// std::nullopt is returned
std::optional<DsCruId> decodeDsCruId(std::string rep);

/// Returns the channel number of a full string representation of (DsCruId,channel)
/// If the string does not correspond to a valid (DsCruId,channel) pair then
/// std::nullopt is returned
//std::optional<uint8_t> decodeChannelId(std::string rep);

std::ostream& operator<<(std::ostream& os, const DsCruId& id);

/// Returns a string representation of the given DsCruId
std::string asString(DsCruId dsId);

/// Extracts the group from the elinkId
std::optional<uint8_t> groupFromElinkId(uint8_t elinkId);

/// Extracts the index from the elinkId
std::optional<uint8_t> indexFromElinkId(uint8_t elinkId);

} // namespace o2::mch::raw

#endif
