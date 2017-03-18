// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Dennis Goldfarb $
// $Authors: Dennis Goldfarb $
// --------------------------------------------------------------------------
//

///////////////////////////

// This one is going to be tested.
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>

///////////////////////////

// More headers

#include <iostream>
#include <iterator>
#include <utility>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IsotopeSplineDB, "$Id$")

/////////////////////////////////////////////////////////////

IsotopeSplineDB* db = 0;
IsotopeSplineDB* db_nullPointer = 0;

START_SECTION(static const IsotopeSplineDB* getInstance())
  db = IsotopeSplineDB::getInstance();
  TEST_NOT_EQUAL(db, db_nullPointer)
END_SECTION

START_SECTION(virtual ~IsotopeSplineDB())
  NOT_TESTABLE
END_SECTION

START_SECTION(UInt inModelBounds(double average_weight, UInt max_depth) const)
  // These are regression tests that will likely break with newly trained splines
  // Just below the minimum mass
  TEST_EQUAL(FALSE, db->inModelBounds(75, 1))
  // Just above the minimum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(76, 1))
  // Just below the maximum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(9992, 1))
  // Just above the maximum mass
  TEST_EQUAL(FALSE, db->inModelBounds(9993, 1))
  // Median mass
  TEST_EQUAL(TRUE,  db->inModelBounds(5000, 1))
  // max_depth can't be 0
  TEST_EQUAL(FALSE, db->inModelBounds(5000, 0))
  // Maximum depth supported by current models
  TEST_EQUAL(TRUE,  db->inModelBounds(5000, 21))
  // Just above the maximum depth
  TEST_EQUAL(FALSE, db->inModelBounds(5000, 22))
  // Both max_depth and mass are out of bounds
  TEST_EQUAL(FALSE, db->inModelBounds(10000, 22))
END_SECTION

START_SECTION(UInt inModelBounds(double average_weight, UInt max_depth, UInt S) const)
  // These are regression tests that will likely break with newly trained splines
  // Just below the minimum mass
  TEST_EQUAL(FALSE, db->inModelBounds(75, 1, 0))
  // Just above the minimum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(76, 1, 0))
  // Just below the maximum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(9999, 1, 0))
  // Just above the maximum mass
  TEST_EQUAL(FALSE, db->inModelBounds(10000, 1, 0))
  // Median mass
  TEST_EQUAL(TRUE,  db->inModelBounds(5000, 1, 0))
  // max_depth can't be 0
  TEST_EQUAL(FALSE, db->inModelBounds(5000, 0, 0))
  // Maximum depth supported by current models
  TEST_EQUAL(TRUE,  db->inModelBounds(5000, 21, 0))
  // Just above the maximum depth
  TEST_EQUAL(FALSE, db->inModelBounds(5000, 22, 0))
  // Both max_depth and mass are out of bounds
  TEST_EQUAL(FALSE, db->inModelBounds(10000, 22, 0))
  // With 10 sulfurs, the minimum mass is higher
  TEST_EQUAL(FALSE,  db->inModelBounds(1000, 1, 10))
  // Just above the minimum mass for a peptide containing 10 sulfurs
  TEST_EQUAL(TRUE,  db->inModelBounds(1049, 1, 10))
  // Just above the maximum mass
  TEST_EQUAL(FALSE,  db->inModelBounds(10000, 1, 10))
  // Just above the maximum number of sulfurs we have splines for
  TEST_EQUAL(FALSE,  db->inModelBounds(5000, 1, 11))
END_SECTION

START_SECTION(IsotopeDistribution estimateFromPeptideWeight(double average_weight, UInt max_depth) const)
  // Regression tests
  IsotopeDistribution iso = db->estimateFromPeptideWeight(100.0, 3);
  TEST_REAL_SIMILAR(iso.begin()->second, 0.94639500435985)

  iso = db->estimateFromPeptideWeight(1000, 3);
  TEST_REAL_SIMILAR(iso.begin()->second, 0.583321621735243)

  iso = db->estimateFromPeptideWeight(9900.0, 3);
  TEST_REAL_SIMILAR(iso.begin()->second, 0.0484190716597988)

  // Comparisons to averagine FFT method
  // max_depth = 0 so FFT method will calculate all isotopes and 0 is out of bounds for the splines
  // so it should default to the averagine method. Results should therefore be identical.
  IsotopeDistribution iso2;
  iso2.estimateFromPeptideWeight(100.0);
  iso = db->estimateFromPeptideWeight(100.0, 0);
  TEST_EQUAL(iso.size(), iso2.size());

  IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso2.begin());
  for (; it1 != iso.end(); ++it1, ++it2)
  {
    TEST_EQUAL(it1->first, it2->first)
    TEST_REAL_SIMILAR(it2->second, it2->second)
  }

  // max_depth = 3 which is in the bounds of the splines. The results should slightly disagree
  IsotopeDistribution iso3(3);
  iso3.estimateFromPeptideWeight(100.0);
  iso = db->estimateFromPeptideWeight(100.0, 3);
  TEST_EQUAL(iso.size(), iso3.size());
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso3.getContainer()[0].second, -0.00334015646431129)


END_SECTION

START_SECTION(IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, UInt S, UInt max_depth) const)
  // Regression tests
  IsotopeDistribution iso, iso2(3);
  // With 0 sulfurs, it should be very unlikely for this tiny peptide to be M+2.
  iso = db->estimateFromPeptideWeightAndS(300.0, 0, 3);
  TEST_REAL_SIMILAR(iso.rbegin()->second, 0.0179838988352706)
  // Results with averagine FFT should be very similar
  iso2.estimateFromPeptideWeightAndS(300.0, 0);
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00152827793431554)
  // With one sulfur, it's more likely that the precursor is M+2 compared to having 0 sulfurs.
  iso = db->estimateFromPeptideWeightAndS(300.0, 1, 3);
  TEST_REAL_SIMILAR(iso.rbegin()->second, 0.0528088116350473)
  // Results with averagine FFT should be very similar
  iso2.estimateFromPeptideWeightAndS(300.0, 1);
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00342043471971454)
  // With two sulfurs, the M+2 isotope is more likely
  iso = db->estimateFromPeptideWeightAndS(300.0, 2, 3);
  TEST_REAL_SIMILAR(iso.rbegin()->second, 0.0858685470961625)
  // Results with averagine FFT should be very similar
  iso2.estimateFromPeptideWeightAndS(300.0, 2);
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00684276895174518)
  // With three sulfurs, the M+2 isotope is even more likely
  iso = db->estimateFromPeptideWeightAndS(300.0, 3, 3);
  TEST_REAL_SIMILAR(iso.rbegin()->second, 0.116873174480882)
  // We're out of bounds of the spline so the averagine FFT method should be identical
  iso2.estimateFromPeptideWeightAndS(300.0, 3);
  TEST_EQUAL(iso.getContainer()[0].second, iso2.getContainer()[0].second)
END_SECTION

START_SECTION(IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor,
                      double average_weight_fragment, const std::set<UInt>& precursor_isotopes) const)
  IsotopeDistribution iso, iso2;
  std::set<UInt> precursor_isotopes;
  // We're isolating the M0 and M+1 precursor isotopes
  precursor_isotopes.insert(0);
  precursor_isotopes.insert(1);
  // These are regression tests, but the results also follow an expected pattern.
  // The differences between the spline and the averagine FFT should be small.

  iso = db->estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second, 0.962439218732925)
  iso2.estimateForFragmentFromPeptideWeight(200.0, 100.0, precursor_isotopes);
  iso2.renormalize();
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00778441741284197)

  iso = db->estimateForFragmentFromPeptideWeight(2000.0, 100.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second, 0.980242098943341)
  iso2.estimateForFragmentFromPeptideWeight(2000.0, 100.0, precursor_isotopes);
  iso2.renormalize();
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00425723273112488)

  iso = db->estimateForFragmentFromPeptideWeight(20000.0, 100.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second, 0.996563468222958)
  iso2.estimateForFragmentFromPeptideWeight(20000.0, 100.0, precursor_isotopes);
  iso2.renormalize();
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.000779946871177439)

  iso = db->estimateForFragmentFromPeptideWeight(2000.0, 1000.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second, 0.7427331415753)
  iso2.estimateForFragmentFromPeptideWeight(2000.0, 1000.0, precursor_isotopes);
  iso2.renormalize();
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.00144216393601726)

  iso = db->estimateForFragmentFromPeptideWeight(20000.0, 1000.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second, 0.955168708458127)
  iso2.estimateForFragmentFromPeptideWeight(20000.0, 1000.0, precursor_isotopes);
  iso2.renormalize();
  TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.000497158581316604)

  // These masses are out of bounds of the splines so the results should be the same as averagine FFT
  iso = db->estimateForFragmentFromPeptideWeight(20000.0, 10000.0, precursor_isotopes);
  iso.renormalize();
  TEST_REAL_SIMILAR(iso.begin()->second,  0.542260764523188)
  iso2.estimateForFragmentFromPeptideWeight(20000.0, 10000.0, precursor_isotopes);
  iso2.renormalize();
  TEST_EQUAL(iso.getContainer()[0].second, iso2.getContainer()[0].second)

  // If the fragment is identical to the precursor, then the distribution
  // should be the same as if it was just a precursor that wasn't isolated.
  iso = db->estimateForFragmentFromPeptideWeight(200.0, 200.0, precursor_isotopes);
  IsotopeDistribution iso_precursor = db->estimateFromPeptideWeight(200.0, 2);
  IsotopeDistribution::ConstIterator it1(iso.begin()), it2(iso_precursor.begin());

  for (; it1 != iso.end(); ++it1, ++it2)
  {
    TEST_EQUAL(it1->first, it2->first)
    TEST_REAL_SIMILAR(it2->second, it2->second)
  }
END_SECTION

START_SECTION(IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor,
                      UInt S_precursor, double average_weight_fragment, UInt S_fragment,
                      const std::set<UInt>& precursor_isotopes) const)
        IsotopeDistribution iso;
        IsotopeDistribution iso2;
        std::set<UInt> precursor_isotopes;
        // We're isolating the M+2 precursor isotopes
        precursor_isotopes.insert(2);
        // These are regression tests, but the results also follow an expected pattern.

        // With 0 sulfurs, it should be somewhat unlikely for the fragment to be M+2.
        iso = db->estimateForFragmentFromPeptideWeightAndS(200.0, 0, 100.0, 0, precursor_isotopes);
        iso.renormalize();
        TEST_REAL_SIMILAR(iso.rbegin()->second, 0.42480315566124)
        iso2.estimateForFragmentFromPeptideWeightAndS(200.0, 0, 100.0, 0, precursor_isotopes);
        iso2.renormalize();
        TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.0693575965376885)
        // With the only sulfur being in the fragment, it's much more likely that the fragment
        // is M+2.
        iso = db->estimateForFragmentFromPeptideWeightAndS(200.0, 1, 100.0, 1, precursor_isotopes);
        iso.renormalize();
        TEST_REAL_SIMILAR(iso.rbegin()->second, 0.867991205760886)
        iso2.estimateForFragmentFromPeptideWeightAndS(200.0, 1, 100.0, 1, precursor_isotopes);
        iso2.renormalize();
        TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.0394515422998836)
        // Both sulfurs are in the fragment, so it's even more likely for the fragment to be M+2
        iso = db->estimateForFragmentFromPeptideWeightAndS(200.0, 2, 100.0, 2, precursor_isotopes);
        iso.renormalize();
        TEST_REAL_SIMILAR(iso.rbegin()->second, 0.928018667348661)
        iso2.estimateForFragmentFromPeptideWeightAndS(200.0, 2, 100.0, 2, precursor_isotopes);
        iso2.renormalize();
        TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.0227855489435035)
        // All 3 sulfurs are in the fragment
        iso = db->estimateForFragmentFromPeptideWeightAndS(200.0, 3, 100.0, 3, precursor_isotopes);
        iso.renormalize();
        TEST_REAL_SIMILAR(iso.rbegin()->second, 0.954837464667933)
        iso2.estimateForFragmentFromPeptideWeightAndS(200.0, 3, 100.0, 3, precursor_isotopes);
        iso2.renormalize();
        TEST_REAL_SIMILAR(iso.getContainer()[0].second - iso2.getContainer()[0].second, 0.0158367426586918)
        // In this case, both the fragment and complementary fragment are out of bounds for the splines
        // so the averagine FFT method should give identical results.
        iso = db->estimateForFragmentFromPeptideWeightAndS(200.0, 3, 150.0, 3, precursor_isotopes);
        iso.renormalize();
        TEST_REAL_SIMILAR(iso.rbegin()->second, 0.975104783221024)
        iso2.estimateForFragmentFromPeptideWeightAndS(200.0, 3, 150.0, 3, precursor_isotopes);
        iso2.renormalize();
        TEST_EQUAL(iso.getContainer()[0].second, iso2.getContainer()[0].second)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
