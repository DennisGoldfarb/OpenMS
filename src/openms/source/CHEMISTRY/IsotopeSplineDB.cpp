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

#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/IsotopeSplineXMLFile.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
    IsotopeSplineDB::IsotopeSplineDB()
    {
      readSplinesFromFile_("CHEMISTRY/IsotopeSplines.xml");
      min_mass_ = 334;
      max_mass_ = 10000;
      max_isotope_ = 20;
    }

    IsotopeSplineDB::~IsotopeSplineDB()
    {
      clear_();
    }

    void IsotopeSplineDB::readSplinesFromFile_(const String& filename)
    {
      String file = File::find(filename);

      IsotopeSplineXMLFile splineFile;
      splineFile.load(file, &this->models);
    }

    void IsotopeSplineDB::approximateIsotopeDistribution(IsotopeDistribution::ContainerType& result, double average_weight, UInt max_isotope)
    {
      result.resize(max_isotope);

      for (UInt isotope = 0; isotope < max_isotope; ++isotope)
      {
        ModelAttributes att(0, isotope);
        double probability = models[att]->eval(average_weight);
        //double probability = std::max(0.0, models[att]->eval(average_weight));
        result[isotope] = make_pair(Size(average_weight + isotope), probability);
      }

    }

    bool IsotopeSplineDB::inModelBounds(double average_weight, UInt max_isotope)
    {
      // Check if masses are in bounds
      if (average_weight < min_mass_ || average_weight >= max_mass_)
      {
        return false;
      }

      // Check if max isotope is in bounds
      if (max_isotope > max_isotope_)
      {
        return false;
      }

      // All checks passed
      return true;
    }

    void IsotopeSplineDB::clear_()
    {
      clearModels_();
    }

    void IsotopeSplineDB::clearModels_()
    {
      for (Iterator itr = models.begin(); itr != models.end(); ++itr)
      {
        delete itr->second;
      }
      models.clear();
    }

    bool operator<(const ModelAttributes &lhs, const ModelAttributes &rhs)
    {
      if (lhs.isotope != rhs.isotope)
      {
        return lhs.isotope < rhs.isotope;
      }
      return lhs.num_sulfur < rhs.num_sulfur;
    }
}
