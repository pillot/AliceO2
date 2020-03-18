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
#include "RawBufferDecoder.h"

// Dans ce code, on récupère un infut aui est un message avec le buffer, on fait tourner le code de base decodeBuffer qui est dans Handlers, et on renvoir un message de sortie (inspiré de FileReader de Andrea)

using namespace o2;
using namespace o2::framework;

class DataDecoderTask
{
  RawBufferDecoder<SampleMode, RDHv4> decoder;

public:
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    // Rien à initialiser
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    bool verbose = false;
    // get the input buffer

    auto msgIn = pc.inputs().get<gsl::span<char>>("readout");
    auto bufferPtrIn = msgIn.data();
    auto sizeIn = msgIn.size();

    if(verbose) std::cout << "We got the input buffer" << std::endl;

    //std::vector<uint8_t> buffer(sizeIn);
    std::vector<uint8_t> buffer((uint8_t*)bufferPtrIn, ((uint8_t*)bufferPtrIn)+sizeIn);

    // Decode the buffer
    size_t outsize;
    char* outbuffer = decoder.decodeBuffer(buffer, outsize);



    /// send the output buffer via DPL

    const int OUT_SIZE = outsize;

    // create the output message
    auto freefct = [](void* data, void* /*hint*/) { free(data); };
    pc.outputs().adoptChunk(Output{ "MCH", "DIGITS", 0 }, (char*)outbuffer, OUT_SIZE, freefct, nullptr);
    //exit(0);
  }

private:
  std::ifstream mInputFile{}; ///< input file
  bool mPrint = false;        ///< print digits

};

// clang-format off
WorkflowSpec defineDataProcessing(const ConfigContext&)
{
  WorkflowSpec specs;

  DataProcessorSpec producer{
    "DataDecoder",
    //Inputs{InputSpec{"readout", "ROUT", "RAWDATA", Lifetime::Timeframe}},
    o2::framework::select("readout:ROUT/RAWDATA"),
    Outputs{OutputSpec{"MCH", "DIGITS", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<DataDecoderTask>()},
    Options{}

  };
  specs.push_back(producer);

  return specs;
}
// clang-format on
