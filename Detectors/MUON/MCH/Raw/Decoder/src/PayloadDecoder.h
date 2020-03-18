// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_RAW_PAYLOAD_DECODER_H
#define O2_MCH_RAW_PAYLOAD_DECODER_H

#include "BareGBTDecoder.h"
#include "DumpBuffer.h"
#include "Headers/RAWDataHeader.h"
#include "MCHRawDecoder/Decoder.h"
#include "MCHRawElecMap/CruLinkId.h"
#include "MakeArray.h"
#include "PayloadDecoder.h"
#include "UserLogicGBTDecoder.h"
#include <cstdlib>
#include <fmt/format.h>
#include <gsl/span>
#include <iostream>

namespace
{
bool hasOrbitJump(uint32_t orb1, uint32_t orb2)
{
  return std::abs(static_cast<long int>(orb1 - orb2)) > 1;
}
} // namespace

namespace o2
{
namespace mch
{
namespace raw
{
/// @brief Decoder for MCH  Raw Data Format.

template <typename RDH, typename GBTDECODER>
class PayloadDecoder
{
 public:
  /// Constructs a decoder
  /// \param channelHandler the handler that will be called for each
  /// piece of sampa data (a SampaCluster, i.e. a part of a time window)
  PayloadDecoder(SampaChannelHandler channelHandler);

  /// decode the buffer
  /// \return the number of bytes used from the buffer
  size_t process(const RDH& rdh, gsl::span<uint8_t> buffer);

  void reset();

 private:
  std::optional<GBTDECODER> mDecoder; //< helper decoders
  SampaChannelHandler mChannelHandler;
  uint32_t mOrbit{0};
};

template <typename RDH, typename GBTDECODER>
PayloadDecoder<RDH, GBTDECODER>::PayloadDecoder(SampaChannelHandler channelHandler)
  : mChannelHandler(channelHandler)
{
}

template <typename RDH, typename GBTDECODER>
size_t PayloadDecoder<RDH, GBTDECODER>::process(const RDH& rdh, gsl::span<uint8_t> buffer)
{
  if (hasOrbitJump(rdhOrbit(rdh), mOrbit)) {
    //std::cout << "Has orbit jump" << std::endl;
    reset();
  } else if (rdhOrbit(rdh) != mOrbit) {
    //std::cout << "diff mOrbit" << std::endl;
  }
  mOrbit = rdhOrbit(rdh);

  int cruId = rdh.feeId & 0xFF;
  int linkId = rdh.linkID;
  uint32_t cruLinkId = o2::mch::raw::encode(o2::mch::raw::CruLinkId(cruId, linkId));
  //std::cout << "[PayloadDecoder] solarId = rdh.feeId = " << int(cruId) << std::endl;
  if(!mDecoder.has_value()) {
    mDecoder.emplace(GBTDECODER(cruLinkId, mChannelHandler));
  }
  return mDecoder.value().append(buffer);
}

template <typename RDH, typename GBTDECODER>
void PayloadDecoder<RDH, GBTDECODER>::reset()
{
  if(mDecoder.has_value()) mDecoder.value().reset();
}

} // namespace raw
} // namespace mch
} // namespace o2

#endif
