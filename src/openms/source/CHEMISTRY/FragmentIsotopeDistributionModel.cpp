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

#include <OpenMS/CHEMISTRY/FragmentIsotopeDistributionModel.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/IsotopeSplineXMLFile.h>

#include <cmath>
#include <algorithm>
#include <iostream>



using namespace std;

namespace OpenMS
{
  FragmentIsotopeDistributionModel::FragmentIsotopeDistributionModel()
  {
    readSplineModelsFromFile_("CHEMISTRY/IsotopeSplines.xml");
    min_fragment_mass = 300;
    min_precursor_mass = 300;
    max_fragment_mass = 8500;
    max_precursor_mass = 8500;
  }

  FragmentIsotopeDistributionModel::~FragmentIsotopeDistributionModel()
  {
    clear_();
  }

  void FragmentIsotopeDistributionModel::readSplineModelsFromFile_(const String& filename) {
    String file = File::find(filename);

    IsotopeSplineXMLFile splineFile;
    splineFile.load(file, &this->models);
  }

  void FragmentIsotopeDistributionModel::approximateIsotopeDistribution(IsotopeDistribution::ContainerType& result, double average_weight_precursor, double average_weight_fragment, const std::vector<UInt>& precursor_isotopes)
  {
    UInt max_depth = *std::max_element(precursor_isotopes.begin(), precursor_isotopes.end())+1;
    result.resize(max_depth);

    for (UInt fragment_isotope = 0; fragment_isotope < max_depth; ++fragment_isotope)
    {
      result[fragment_isotope] = make_pair(Size(average_weight_fragment + fragment_isotope), 0);
    }

    for (std::vector<UInt>::const_iterator precursor_isotope_itr = precursor_isotopes.begin(); precursor_isotope_itr != precursor_isotopes.end(); ++precursor_isotope_itr)
    {
      for (UInt fragment_isotope = 0; fragment_isotope <= *precursor_isotope_itr; ++fragment_isotope)
      {
          // get model index
          //UInt model_index = getModelIndex(average_weight_precursor, average_weight_fragment, *precursor_isotope_itr, fragment_isotope);
          ModelAttributes att(0, 0, 0, 0, *precursor_isotope_itr, fragment_isotope);

          // add contribution for this fragment isotope from this precursor isotope
          result[fragment_isotope].second += models[att]->evaluate_model(average_weight_precursor, average_weight_fragment);
      }
    }
  }

  bool FragmentIsotopeDistributionModel::inModelBounds(double average_weight_precursor, double average_weight_fragment, const std::vector<UInt>& precursor_isotopes)
  {
    // Check if masses are in bounds
    if (average_weight_precursor < min_precursor_mass || average_weight_precursor > max_precursor_mass
          || average_weight_fragment < min_fragment_mass || average_weight_fragment > max_fragment_mass)
    {
      return false;
    }

    // Check if isolated isotopes are in bounds
    for (std::vector<UInt>::const_iterator itr = precursor_isotopes.begin(); itr != precursor_isotopes.end(); ++itr)
    {
      if (*itr > max_isotope)
      {
        return false;
      }
    }

    // All checks passed
    return true;
  }

  UInt FragmentIsotopeDistributionModel::getModelIndex(double precursor_mass, double fragment_mass, UInt precursor_isotope, UInt fragment_isotope)
  {
    UInt models_per_precursor_isotope = max_isotope + 1;
    UInt models_per_fragment_mass = (max_isotope + 1) * models_per_precursor_isotope;
    UInt models_per_precursor_mass = ((max_fragment_mass - min_fragment_mass) / fragment_mass_step) * models_per_fragment_mass;

    UInt precursor_mass_index = floor((precursor_mass-min_precursor_mass) / precursor_mass_step);
    UInt fragment_mass_index = floor((fragment_mass-min_fragment_mass) / fragment_mass_step);

    return (precursor_mass_index * models_per_precursor_mass)
           + (fragment_mass_index * models_per_fragment_mass)
           + (precursor_isotope * models_per_precursor_isotope)
           + fragment_isotope;
  }

  void FragmentIsotopeDistributionModel::clear_() {
    clearModels_();
  }

  void FragmentIsotopeDistributionModel::clearModels_() {
    for (Iterator itr = models.begin(); itr != models.end(); ) {
      itr = models.erase(itr);
    }
    models.clear();
  }

}
