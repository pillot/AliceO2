// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "PageParser.h"
#include "BareGBTDecoder.h"
#include "UserLogicGBTDecoder.h"
#include "Headers/RAWDataHeader.h"
#include "MCHRawCommon/DataFormats.h"

namespace o2::mch::raw
{

  template <typename FORMAT, typename CHARGESUM>
  struct GBTDecoderTrait {
    using type = void;
  };

  template <typename CHARGESUM>
  struct GBTDecoderTrait<BareFormat, CHARGESUM> {
    using type = BareGBTDecoder<CHARGESUM>;
  };

  template <typename CHARGESUM>
  struct GBTDecoderTrait<UserLogicFormat, CHARGESUM> {
    using type = UserLogicGBTDecoder<CHARGESUM>;
  };

template <typename CHARGESUM, typename RDH>
Decoder createDecoder(RawDataHeaderHandler<RDH> rdhHandler, SampaChannelHandler channelHandler)
{
  using BGBTDecoder = typename GBTDecoderTrait<BareFormat, CHARGESUM>::type;
  using BPAYLOADDECODER = PayloadDecoder<RDH, BGBTDecoder>;
  using ULGBTDecoder = typename GBTDecoderTrait<UserLogicFormat, CHARGESUM>::type;
  using ULPAYLOADDECODER = PayloadDecoder<RDH, ULGBTDecoder>;
  return [rdhHandler, channelHandler](gsl::span<uint8_t> buffer) -> DecoderStat {
    static PageParser<RDH, BPAYLOADDECODER, ULPAYLOADDECODER> mPageParser(rdhHandler, BPAYLOADDECODER(channelHandler), ULPAYLOADDECODER(channelHandler));
    return mPageParser.parse(buffer);
  };
}

std::ostream& operator<<(std::ostream& out, const DecoderStat& decStat)
{
  out << fmt::format("Nof orbits seen {} - Nof orbits jumps {}",
                     decStat.nofOrbitSeen,
                     decStat.nofOrbitJumps);
  return out;
}

// define only the specialization we use

using RDHv4 = o2::header::RAWDataHeaderV4;

template Decoder createDecoder<SampleMode, RDHv4>(RawDataHeaderHandler<RDHv4>, SampaChannelHandler);
template Decoder createDecoder<ChargeSumMode, RDHv4>(RawDataHeaderHandler<RDHv4>, SampaChannelHandler);
} // namespace o2::mch::raw
