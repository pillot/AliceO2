// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file MathiesonOriginal.cxx
/// \brief Original implementation of the Mathieson function
///
/// \author Philippe Pillot, Subatech

#include "MCHBase/MathiesonOriginal.h"

#include <TMath.h>

namespace o2
{
namespace mch
{

//_________________________________________________________________________________________________
void MathiesonOriginal::setSqrtKx3AndDeriveKx2Kx4(double sqrtKx3)
{
  /// set the Mathieson parameter sqrt(K3) in x direction, perpendicular to the wires,
  /// and derive the Mathieson parameters K2 and K4 in the same direction
  mSqrtKx3 = sqrtKx3;
  mKx2 = TMath::Pi() / 2. * (1. - 0.5 * mSqrtKx3);
  mKx4 = 0.25 / TMath::ATan(mSqrtKx3);
}

//_________________________________________________________________________________________________
void MathiesonOriginal::setSqrtKy3AndDeriveKy2Ky4(double sqrtKy3)
{
  /// set the Mathieson parameter sqrt(K3) in y direction, along the wires,
  /// and derive the Mathieson parameters K2 and K4 in the same direction
  mSqrtKy3 = sqrtKy3;
  mKy2 = TMath::Pi() / 2. * (1. - 0.5 * mSqrtKy3);
  mKy4 = 0.25 / TMath::ATan(mSqrtKy3);
}

//_________________________________________________________________________________________________
double MathiesonOriginal::integrate(double xMin, double yMin, double xMax, double yMax) const
{
  /// integrate the Mathieson over x and y in the given area

  xMin *= mInversePitch;
  xMax *= mInversePitch;
  yMin *= mInversePitch;
  yMax *= mInversePitch;
  //
  // The Mathieson function
  double uxMin = mSqrtKx3 * TMath::TanH(mKx2 * xMin);
  double uxMax = mSqrtKx3 * TMath::TanH(mKx2 * xMax);

  double uyMin = mSqrtKy3 * TMath::TanH(mKy2 * yMin);
  double uyMax = mSqrtKy3 * TMath::TanH(mKy2 * yMax);

  return 4. * mKx4 * (TMath::ATan(uxMax) - TMath::ATan(uxMin)) *
         mKy4 * (TMath::ATan(uyMax) - TMath::ATan(uyMin));
}

} // namespace mch
} // namespace o2
