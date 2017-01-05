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

      max_isotope_ = this->models.size()-1;
    }

    void IsotopeSplineDB::approximateIsotopeDistribution(IsotopeDistribution::ContainerType& result, double average_weight, UInt max_isotope)
    {
      UInt max_depth = max_isotope;
      result.resize(max_depth);

      for (UInt isotope = 0; isotope < max_depth; ++isotope)
      {
        double probability = models[isotope].eval(average_weight);
        //double probability = std::max(0.0, models[att]->eval(average_weight));
        result[isotope] = make_pair(Size(average_weight + isotope), probability);
      }

    }

    bool IsotopeSplineDB::inModelBounds(double average_weight, UInt max_isotope)
    {
      // Check if max isotope is in bounds
      if (max_isotope > max_isotope_)
      {
        return false;
      }

      // Check if masses are in bounds
      for (UInt isotope = 0; isotope < max_isotope; ++isotope)
      {
        if (!models[isotope].inBounds(average_weight))
        {
          return false;
        }
      }

      // All checks passed
      return true;
    }

    IsotopeDistribution IsotopeSplineDB::estimateFromPeptideWeight(double average_weight, UInt max_depth)
    {
      IsotopeDistribution id(max_depth);
      // Check if the splines can completely handle this request
      if (getInstance()->inModelBounds(average_weight, getInstance()->max_isotope_))
      {
        IsotopeDistribution::ContainerType result;
        IsotopeSplineDB::getInstance()->approximateIsotopeDistribution(result, average_weight, max_isotope_);
        id.set(result);
      }
      else
      {
        id.estimateFromPeptideWeight(average_weight);
      }

      return id;
    }

    IsotopeDistribution IsotopeSplineDB::estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, const std::set<UInt>& precursor_isotopes)
    {
      UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end()) + 1;

      IsotopeDistribution id_fragment(max_depth);
      id_fragment.estimateFromPeptideWeightFast(average_weight_fragment);

      IsotopeDistribution id_comp_fragment(max_depth);
      id_comp_fragment.estimateFromPeptideWeightFast(average_weight_precursor - average_weight_fragment);

      IsotopeDistribution result(max_depth);
      result.calcFragmentIsotopeDist(id_fragment, id_comp_fragment, precursor_isotopes);

      return result;
    }


      void IsotopeSplineDB::clear_()
    {
      clearModels_();
    }

    void IsotopeSplineDB::clearModels_()
    {
      models.clear();
    }
}
