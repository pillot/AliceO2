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
/// \file    runFileReader.cxx
/// \author  Andrea Ferrero
///
/// \brief This is an executable that reads a data file from disk and sends the data to QC via DPL.
///
/// This is an executable that reads a data file from disk and sends the data to QC via the Data Processing Layer.
/// It can be used as a data source for QC development. For example, one can do:
/// \code{.sh}
/// o2-qc-run-file-reader --infile=some_data_file | o2-qc --config json://${QUALITYCONTROL_ROOT}/etc/your_config.json
/// \endcode
///

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "MCHBase/Digit.h"
#include "DataDecoderSpec.h"
#include "Handlers.h"

// Dans ce code, on récupère un infut aui est un message avec le buffer, on fait tourner le code de base decodeBuffer qui est dans Handlers, et on renvoir un message de sortie (inspiré de FileReader de Andrea)

using namespace o2;
using namespace o2::framework;

class DigitReaderTask
{
 public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
// Rien à initialiser
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {

        // get the input buffer
      
      for (auto&& input : pc.inputs()) {
         std::cout << "run RawDataProcessor: input " << input.spec->binding << std::endl;
      }
      
      std::cout << "A" << std::endl;
      
        auto msgIn = pc.inputs().get<gsl::span<char>>("readout");
        auto bufferPtrIn = msgIn.data();
        auto sizeIn = msgIn.size();
      
      std::cout << "We got the input buffer" << std::endl;
      
      std::vector<uint8_t> buffer(sizeIn);
      
      
        // get header info and check message consistency
//                    if (sizeIn < SSizeOfDigitBlock) {
//                      throw out_of_range("missing DigitBlock");
//                    }
//                    auto digitBlock(reinterpret_cast<const DigitBlock*>(bufferPtrIn));
//                    bufferPtrIn += SSizeOfDigitBlock;
//                    sizeIn -= SSizeOfDigitBlock;
//                    if (digitBlock->header.fRecordWidth != SSizeOfDigitStruct) {
//                      throw length_error("incorrect size of digits. Corrupted message?");
//                    }
//                    if (sizeIn != digitBlock->header.fNrecords * SSizeOfDigitStruct) {
//                      throw length_error("incorrect payload");
//                    }

        // load the digits to get the fired pads
//                    auto digits(reinterpret_cast<const DigitStruct*>(bufferPtrIn));

        // Add digits conversion here
//                    mPreClusterFinder.loadDigits(digits, digitBlock->header.fNrecords);

        // Decode the buffer
      size_t outsize;
        char* outbuffer = decodeBuffer<BareFormat, ChargeSumMode, RDHv4>(buffer, outsize);

      

    /// send the output buffer via DPL

    const int OUT_SIZE = outsize;

    // create the output message
    auto freefct = [](void* data, void* /*hint*/) { free(data); };
    pc.outputs().adoptChunk(Output{ "MCH", "DIGITS" }, (char*)outbuffer, OUT_SIZE, freefct, nullptr);
  }

 private:
  std::ifstream mInputFile{}; ///< input file
  bool mPrint = false;        ///< print digits
  
};

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getDigitReaderSpec()
{
  return DataProcessorSpec{
    "DigitReader",
    Inputs{InputSpec{"readout", "ROUT", "RAWDATA", Lifetime::Timeframe}},
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{ adaptFromTask<DigitReaderTask>() },
    Options{ { "infile", VariantType::String, "data.raw", { "input file name" } } }
  };
}

// clang-format off
WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  // The producer to generate some data in the workflow
    
    // Outputs{ { { "digits" }, { "MCH", "DIGITS" } } },
    
  DataProcessorSpec producer{
    "DigitReader",
    Inputs{InputSpec{"readout", "ROUT", "RAWDATA", Lifetime::Timeframe}},
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<DigitReaderTask>()},
    Options{ { "infile", VariantType::String, "", { "input digits" } } }
      
  };
  specs.push_back(producer);

  return specs;
}
// clang-format on
