#ifndef O2_MCH_DATADECODERSPEC_H_
#define O2_MCH_DATADECODERSPEC_H_

#include <random>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "Framework/CallbackService.h"
#include "Framework/ControlService.h"
#include "Framework/Task.h"
#include "Framework/runDataProcessing.h"
#include "MCHBase/Digit.h"

using namespace o2;
using namespace o2::framework;


o2::framework::DataProcessorSpec getDigitReaderSpec();
WorkflowSpec defineDataProcessing(const ConfigContext&);

#endif

