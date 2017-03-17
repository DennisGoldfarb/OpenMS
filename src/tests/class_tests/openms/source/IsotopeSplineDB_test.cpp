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
  TEST_EQUAL(FALSE, db->inModelBounds(333, 1))
  // Just above the minimum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(334, 1))
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
  TEST_EQUAL(FALSE, db->inModelBounds(333, 1, 0))
  // Just above the minimum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(334, 1, 0))
  // Just below the maximum mass
  TEST_EQUAL(TRUE,  db->inModelBounds(9992, 1, 0))
  // Just above the maximum mass
  TEST_EQUAL(FALSE, db->inModelBounds(9993, 1, 0))
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

  TEST_EQUAL(TRUE,  db->inModelBounds(1000, 1, 10))
  TEST_EQUAL(TRUE,  db->inModelBounds(9000, 1, 10))
  TEST_EQUAL(TRUE,  db->inModelBounds(9000, 20, 10))

  TEST_EQUAL(FALSE,  db->inModelBounds(1000, 1, 11))
  TEST_EQUAL(FALSE,  db->inModelBounds(9000, 1, 11))
  TEST_EQUAL(FALSE,  db->inModelBounds(9000, 20, 11))

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
