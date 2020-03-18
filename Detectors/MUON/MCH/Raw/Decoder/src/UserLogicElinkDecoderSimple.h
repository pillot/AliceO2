// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_MCH_RAW_USER_LOGIC_ELINK_DECODER_SIMPLE_H
#define O2_MCH_RAW_USER_LOGIC_ELINK_DECODER_SIMPLE_H

#include "Assertions.h"
#include "MCHRawCommon/DataFormats.h"
#include "MCHRawCommon/SampaHeader.h"
#include "MCHRawDecoder/SampaChannelHandler.h"
#include <bitset>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>

namespace o2::mch::raw
{

  enum decode_state_t
  {
    DECODE_STATE_UNKNOWN,
    DECODE_STATE_SYNC_FOUND,
    DECODE_STATE_HEADER_FOUND,
    DECODE_STATE_CSIZE_FOUND,
    DECODE_STATE_CTIME_FOUND,
    DECODE_STATE_SAMPLE_FOUND,
    DECODE_STATE_PACKET_END
  };


/// @brief Main element of the MCH Bare Raw Data Format decoder.
///
/// A UserLogicElinkDecoder manages the bit stream for one Elink.
///
/// Bits coming from parts of the GBT words are added to the Elink using the
/// append() method and each time a SampaCluster is decoded,
/// it is passed to the SampaChannelHandler for further processing (or none).
///
/// \nosubgrouping
///
template <typename CHARGESUM>
class UserLogicElinkDecoder
{
 public:
  /// Constructor.
  /// \param dsId the (electronic) id of the dual sampa this elink
  /// is connected  to
  /// \param sampaChannelHandler a callable that is passed
  /// each SampaCluster that will be decoded
  UserLogicElinkDecoder(DsCruId dsId, SampaChannelHandler sampaChannelHandler);

  /** @name Main interface 
  */
  ///@{

  /// Append two bits (from the same dual sampa, one per sampa) to the Elink.
  void append(uint64_t data);
  ///@}

  // /// linkId is the GBT id this Elink is part of
  // uint8_t linkId() const;

  /** @name Methods for testing
    */
  ///@{

  /// Current number of bits we're holding
  int len() const;

  /// Reset our internal bit stream, and the sync status
  /// i.e. assume the sync has to be found again
  void reset();
  ///@}

 private:
  /// The possible states we can be in
  enum class State : int {
    LookingForSync,    //< we've not found a sync yet
    LookingForHeader,  //< we've looking for a 50-bits header
    ReadingNofSamples, //< we're (about to) read nof of samples
    ReadingTimestamp,  //< we're (about to) read a timestamp (for the current cluster)
    ReadingSample,     //< we're (about to) read a sample (for the current cluster)
    ReadingClusterSum  //< we're (about to) read a chargesum (for the current cluster)
  };

  std::string name(State state) const;
  void append10Bits(uint16_t data);
  void changeState(State newState, int newCheckpoint);
  void changeToReadingData();
  void clear(int checkpoint);
  void findSync();
  void handlReadClusterSum();
  void handleHeader();
  void handleReadClusterSum();
  void handleReadData();
  void handleReadSample();
  void handleReadTimestamp();
  void oneLess10BitWord();
  void process();
  void sendCluster();
  void softReset();

  template <typename T>
  friend std::ostream& operator<<(std::ostream& os, const o2::mch::raw::UserLogicElinkDecoder<T>& e);

 private:
  DsCruId mDsId;
  SampaChannelHandler mSampaChannelHandler; //< The callable that will deal with the SampaCluster objects we decode
  SampaHeader mSampaHeader;                 //< Current SampaHeader
  uint64_t mBitBuffer;                      //< Our internal bit stream buffer
  /** @name internal global counters
    */

  ///@{
  uint64_t mNofSync;               //< Number of SYNC words we've seen so far
  uint64_t mNof10BitSeen;          //< Total number of 10 bits seen
  uint64_t mNofHeaderSeen;         //< Total number of headers seen
  uint64_t mNofHammingErrors;      //< Total number of hamming errors seen
  uint64_t mNofHeaderParityErrors; //< Total number of header parity errors seen
  ///@}

  uint64_t mCheckpoint;           //< mask of the next state transition check to be done in process()
  uint16_t mNof10BitsWordsToRead; //< number of 10 bits words to be read

  uint16_t mNofSamples;
  uint16_t mTimestamp;
  std::vector<uint16_t> mSamples;
  uint32_t mClusterSum;
  uint64_t mMask;

  State mState; //< the state we are in

  bool verbose;
};

constexpr int SAMPA_HEADERSIZE = 50;

//FIXME: probably needs the GBT id as well here ?
template <typename CHARGESUM>
UserLogicElinkDecoder<CHARGESUM>::UserLogicElinkDecoder(DsCruId dsId,
                                              SampaChannelHandler sampaChannelHandler)
  : mDsId{dsId},
    mSampaChannelHandler{sampaChannelHandler},
    mSampaHeader{},
    mBitBuffer{},
    mNofSync{},
    mNof10BitSeen{},
    mNofHeaderSeen{},
    mNofHammingErrors{},
    mNofHeaderParityErrors{},
    mCheckpoint{SAMPA_HEADERSIZE},
    mNof10BitsWordsToRead{},
    mNofSamples{},
    mTimestamp{},
    mSamples{},
    mClusterSum{},
    mState{State::LookingForSync},
    mMask{0},
    verbose(false)
{
      //if(dsId.elinkId() == 2) verbose = true;
      //else verbose = false;
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::append10Bits(uint16_t data)
{
  mNof10BitSeen++;

  mBitBuffer += static_cast<uint64_t>(data) << mMask;
  mMask += 10;

  //std::cout << fmt::format("[append10Bits] data={:08X} bitBuffer={:08X} mask={} checkpoint={}\n", data, mBitBuffer, mMask, mCheckpoint);
  if (mMask == mCheckpoint) {
    process();
  }
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::append(uint64_t data)
{
  int nch = (data >> 53) & 0x7FF;
  int link_id = (data >> 59) & 0x1F;
  int ds_id = (data >> 53) & 0x3F;
  //if(link_id != 1) continue;
  int is_incomplete = (data >> 52) & 0x1;
  int err_code = (data >> 50) & 0x3;

  append10Bits(data&0x3FF);       if(mState == State::LookingForHeader && is_incomplete ) return;
  append10Bits((data>>10)&0x3FF); if(mState == State::LookingForHeader && is_incomplete ) return;
  append10Bits((data>>20)&0x3FF); if(mState == State::LookingForHeader && is_incomplete ) return;
  append10Bits((data>>30)&0x3FF); if(mState == State::LookingForHeader && is_incomplete ) return;
  append10Bits((data>>40)&0x3FF); if(mState == State::LookingForHeader && is_incomplete ) return;
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::changeState(State newState, int newCheckpoint)
{
  mState = newState;
  if(verbose)
    std::cout << fmt::format("[changeState {}-{}] changed to ", mDsId.elinkGroupId(), mDsId.elinkIndexInGroup())
  << name(mState) << fmt::format(", new checkpoint={}\n", newCheckpoint);
  clear(newCheckpoint);
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::clear(int checkpoint)
{
  mBitBuffer = 0;
  mCheckpoint = checkpoint;
  mMask = 0;
}

/// findSync checks if the last 50 bits of the bit stream
/// match the Sampa SYNC word.
///
/// - if they are then reset the bit stream and sets the checkpoint to 50 bits
/// - if they are not then pop the first bit out
template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::findSync()
{
  const uint64_t sync = sampaSync().uint64();
  assert(mState == State::LookingForSync);
  if(verbose) std::cout << fmt::format("[findSync] bitBuffer={:08X} mask={} checkpoint={}\n", mBitBuffer, mMask, mCheckpoint);
  if (mBitBuffer != sync) {
    mBitBuffer >>= 10;
    mMask -= 10;
    return;
  }
  if(verbose) std::cout << fmt::format("[findSync] SYNC found\n");
  changeState(State::LookingForHeader, SAMPA_HEADERSIZE);
  mNofSync++;
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::handleHeader()
{
  assert(mState == State::LookingForHeader);

  if(verbose) std::cout << fmt::format("[handleHeader] bitBuffer={:08X} mask={} checkpoint={}\n", mBitBuffer, mMask, mCheckpoint);

  // check if the current 50-bit word is a SYNC, and skip it if that's the case
  const uint64_t sync = sampaSync().uint64();
  if (mBitBuffer == sync) {
    mNofSync++;
    softReset();
    return;
  }


  mSampaHeader.uint64(mBitBuffer);

  ++mNofHeaderSeen;

  if (mSampaHeader.hasError()) {
    ++mNofHammingErrors;
  }

  switch (mSampaHeader.packetType()) {
    case SampaPacketType::DataTruncated:
    case SampaPacketType::DataTruncatedTriggerTooEarly:
    case SampaPacketType::DataTriggerTooEarly:
    case SampaPacketType::DataTriggerTooEarlyNumWords:
    case SampaPacketType::DataNumWords:
      // data with a problem is still data, i.e. there will
      // probably be some data words to read in...
      // so we fallthrough the simple Data case
    case SampaPacketType::Data:
      mNof10BitsWordsToRead = mSampaHeader.nof10BitWords();
      if(verbose) std::cout << fmt::format("[handleHeader] data header found, nof10BitsWordsToRead={}\n", mNof10BitsWordsToRead);
      changeState(State::ReadingNofSamples, 10);
      break;
    case SampaPacketType::HeartBeat:
      fmt::printf("UserLogicElinkDecoder %d: HEARTBEAT found. Should be doing sth about it ?\n", mDsId);
      softReset();
      break;
    default:
      throw std::logic_error("that should not be possible");
      break;
  }
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::handleReadClusterSum()
{
  mClusterSum = mBitBuffer;
  oneLess10BitWord();
  oneLess10BitWord();
  sendCluster();
  if (mNof10BitsWordsToRead) {
    changeState(State::ReadingNofSamples, 10);
  } else {
    changeState(State::LookingForHeader, SAMPA_HEADERSIZE);
  }
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::handleReadData()
{
  assert(mState == State::ReadingTimestamp || mState == State::ReadingSample);
  if (mState == State::ReadingTimestamp) {
    mTimestamp = mBitBuffer;
  }
  oneLess10BitWord();
  changeToReadingData();
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::handleReadSample()
{
  mSamples.push_back(mBitBuffer);
  if (mNofSamples > 0) {
    --mNofSamples;
  }
  oneLess10BitWord();
  if(verbose) std::cout << fmt::format("[handleReadSample] mNofSamples={} mNof10BitsWordsToRead={} samples.size={}\n",
      mNofSamples, mNof10BitsWordsToRead, mSamples.size());
  if (mNofSamples) {
    handleReadData();
  } else {
    sendCluster();
    if (mNof10BitsWordsToRead) {
      changeState(State::ReadingNofSamples, 10);
    } else {
      changeState(State::LookingForHeader, SAMPA_HEADERSIZE);
    }
  }
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::handleReadTimestamp()
{
  assert(mState == State::ReadingNofSamples);
  oneLess10BitWord();
  mNofSamples = mBitBuffer;
  changeState(State::ReadingTimestamp, 10);
}

template <typename CHARGESUM>
int UserLogicElinkDecoder<CHARGESUM>::len() const
{
  return static_cast<int>(std::floor(log2(1.0 * mMask)) + 1);
}

template <typename CHARGESUM>
std::string UserLogicElinkDecoder<CHARGESUM>::name(State s) const
{
  switch (s) {
    case State::LookingForSync:
      return "LookingForSync";
      break;
    case State::LookingForHeader:
      return "LookingForHeader";
      break;
    case State::ReadingNofSamples:
      return "ReadingNofSamples";
      break;
    case State::ReadingTimestamp:
      return "ReadingTimestamp";
      break;
    case State::ReadingSample:
      return "ReadingSample";
      break;
    case State::ReadingClusterSum:
      return "ReadingClusterSum";
      break;
  };
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::oneLess10BitWord()
{
  if (mNof10BitsWordsToRead > 0) {
    --mNof10BitsWordsToRead;
  }
}

/// process the bit stream content.
template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::process()
{
  if(verbose) std::cout << fmt::format("[process {}-{}] bitBuffer={:08X} state=", mDsId.elinkGroupId(), mDsId.elinkIndexInGroup(), mBitBuffer) << name(mState) << std::endl;
  switch (mState) {
    case State::LookingForSync:
      findSync();
      break;
    case State::LookingForHeader:
      handleHeader();
      break;
    case State::ReadingNofSamples:
      handleReadTimestamp();
      break;
    case State::ReadingTimestamp:
      handleReadData();
      break;
    case State::ReadingSample:
      handleReadSample();
      break;
    case State::ReadingClusterSum:
      handleReadClusterSum();
      break;
  }
};

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::softReset()
{
  clear(SAMPA_HEADERSIZE);
}

template <typename CHARGESUM>
void UserLogicElinkDecoder<CHARGESUM>::reset()
{
  softReset();
  mState = State::LookingForSync;
}

template <typename CHARGESUM>
std::ostream& operator<<(std::ostream& os, const o2::mch::raw::UserLogicElinkDecoder<CHARGESUM>& e)
{
  os << fmt::format("ID{:2d} cruId {:2d} sync {:6d} cp 0x{:6x} mask 0x{:6x} state {:17s} len {:6d} nseen {:6d} errH {:6} errP {:6} head {:6d} n10w {:6d} nsamples {:6d} mode {} bbuf {:s}",
                    e.mLinkId, e.mCruId, e.mNofSync, e.mCheckpoint, e.mMask,
                    e.name(e.mState),
                    e.len(), e.mNofBitSeen,
                    e.mNofHeaderSeen,
                    e.mNofHammingErrors,
                    e.mNofHeaderParityErrors,
                    e.mNof10BitsWordsToRead,
                    e.mNofSamples,
                    (e.mClusterSumMode ? "CLUSUM" : "SAMPLE"),
                    bitBufferString(e.mBitBuffer, e.mMask));
  return os;
}

//uint8_t channelNumber(const SampaHeader& sh)
//{
//  return sh.channelAddress() + (sh.chipAddress() % 2) * 32;
//}

template <>
void UserLogicElinkDecoder<ChargeSumMode>::sendCluster()
{
  if (mSampaChannelHandler) {
    mSampaChannelHandler(mDsId,
                         channelNumber(mSampaHeader),
                         SampaCluster(mTimestamp, mClusterSum));
  }
}

template <>
void UserLogicElinkDecoder<SampleMode>::sendCluster()
{
  if (mSampaChannelHandler) {
    mSampaChannelHandler(mDsId,
                         channelNumber(mSampaHeader),
                         SampaCluster(mTimestamp, mSamples));
  }
  mSamples.clear();
}

template <>
void UserLogicElinkDecoder<ChargeSumMode>::changeToReadingData()
{
  changeState(State::ReadingClusterSum, 20);
}

template <>
void UserLogicElinkDecoder<SampleMode>::changeToReadingData()
{
  changeState(State::ReadingSample, 10);
}

} // namespace o2::mch::raw

#endif
