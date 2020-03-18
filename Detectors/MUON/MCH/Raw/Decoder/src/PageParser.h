// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_RAW_PAGEPARSER_H
#define O2_MCH_RAW_PAGEPARSER_H

#include <gsl/span>
#include "MCHRawDecoder/RawDataHeaderHandler.h"
#include "PayloadDecoder.h"

#define O2_MCH_MAX_CRU_ID 63
#define O2_MCH_MAX_LINK_ID 11

namespace o2::mch::raw
{

template <typename RDH, typename BAREPAYLOADDECODER, typename ULPAYLOADDECODER>
class PageParser
{
 public:
  PageParser(RawDataHeaderHandler<RDH> rdhHandler, BAREPAYLOADDECODER bareDecoder, ULPAYLOADDECODER userLogicDecoder);

  DecoderStat parse(gsl::span<uint8_t> buffer);

 private:
  RawDataHeaderHandler<RDH> mRdhHandler;
  std::map<uint32_t, BAREPAYLOADDECODER> mBareDecoders;
  std::map<uint32_t, ULPAYLOADDECODER>   mUserLogicDecoders;
  BAREPAYLOADDECODER mBareDecoder;
  ULPAYLOADDECODER mUserLogicDecoder;
  DecoderStat mStats;
};

template <typename RDH, typename BAREPAYLOADDECODER, typename ULPAYLOADDECODER>
PageParser<RDH, BAREPAYLOADDECODER, ULPAYLOADDECODER>::PageParser(RawDataHeaderHandler<RDH> rdhHandler, BAREPAYLOADDECODER bareDecoder, ULPAYLOADDECODER userLogicDecoder)
  : mRdhHandler(rdhHandler), mBareDecoder(bareDecoder), mUserLogicDecoder(userLogicDecoder), mStats{}
{
}

template <typename RDH, typename BAREPAYLOADDECODER, typename ULPAYLOADDECODER>
DecoderStat PageParser<RDH, BAREPAYLOADDECODER, ULPAYLOADDECODER>::parse(gsl::span<uint8_t> buffer)
{
  RDH originalRDH;
  const size_t nofRDHWords = sizeof(originalRDH);
  size_t index{0};
  uint64_t nbytes{0};

  while ((index + nofRDHWords) < buffer.size()) {
    originalRDH = createRDH<RDH>(buffer.subspan(index, nofRDHWords));
    if (!isValid(originalRDH)) {
      std::cout << "Got an invalid RDH\n";
      impl::dumpBuffer(buffer.subspan(index, nofRDHWords));
      return mStats;
    }
    auto rdhOpt = mRdhHandler(originalRDH);
    if (!rdhOpt.has_value()) {
      break;
    }
    auto rdh = rdhOpt.value();

    int payloadSize = rdhPayloadSize(rdh);
    size_t n = static_cast<size_t>(payloadSize);
    if (n) {
      size_t pos = static_cast<size_t>(index + nofRDHWords);

      int cruId = rdh.feeId & 0xFF;
      int linkId = rdh.linkID;
      uint32_t cruLinkId = o2::mch::raw::encode(o2::mch::raw::CruLinkId(cruId, linkId));

      if( linkId == 15 ) {
        auto c = mUserLogicDecoders.find(cruLinkId);
        if (c == mUserLogicDecoders.end()) {
          mUserLogicDecoders.emplace(cruLinkId, mUserLogicDecoder);
          c = mUserLogicDecoders.find(cruLinkId);
        }
        c->second.process(rdh, buffer.subspan(pos, n));
      } else {
        auto c = mBareDecoders.find(cruLinkId);
        if (c == mBareDecoders.end()) {
          mBareDecoders.emplace(cruLinkId, mBareDecoder);
          c = mBareDecoders.find(cruLinkId);
        }
        c->second.process(rdh, buffer.subspan(pos, n));
      }
      nbytes += n + nofRDHWords;
    }
    index += rdh.offsetToNext;
  }
  mStats.nofBytesUsed += nbytes;
  return mStats;
}

} // namespace o2::mch::raw

#endif
