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

#include <OpenMS/MATH/MISC/TensorProductSpline.h>
#include <OpenMS/FORMAT/Base64.h>

#include <iostream>


using namespace std;

namespace OpenMS
{

    float TensorProductSpline::evaluate_model(float pmass, float fmass) {
      if (fmass >= pmass) return -1;
      if (fmass <= breaks_fragment_mass[0] || fmass >= breaks_fragment_mass[breaks_fragment_mass.size()-2]) return -1;
      if (pmass <= breaks_precursor_mass[0] || pmass >= breaks_precursor_mass[breaks_precursor_mass.size()-2]) return -1;

      // find index in precursor breaks
      int precursor_index = std::lower_bound(breaks_precursor_mass.begin(), breaks_precursor_mass.end(), pmass)-breaks_precursor_mass.begin()-1;

      // find index in fragment breaks
      int fragment_index = std::lower_bound(breaks_fragment_mass.begin(), breaks_fragment_mass.end(), fmass)-breaks_fragment_mass.begin()-1;

      //if (breaks_fragment_mass[fragment_index] >= breaks_precursor_mass[precursor_index+1]) return -1;

      // do the math
      float* c = coefficients[precursor_index][fragment_index];

      float x = fmass-breaks_fragment_mass[fragment_index];
      float y = pmass-breaks_precursor_mass[precursor_index];

      float v = c[15] + y*(c[14] + y*(c[13] + y*c[12]))
                + x*(c[11] + y*(c[10] + y*(c[9] + y*c[8]))
                     + x*(c[7] + y*(c[6] + y*(c[5] + y*c[4]))
                          + x*(c[3] + y*(c[2] + y*(c[1] + y*c[0])) )));

      return v;
    }

    bool operator<(const ModelAttributes &lhs, const ModelAttributes &rhs) {
      if (lhs.precursor_isotope != rhs.precursor_isotope) return lhs.precursor_isotope < rhs.precursor_isotope;
      if (lhs.fragment_isotope != rhs.fragment_isotope) return lhs.fragment_isotope < rhs.fragment_isotope;
      if (lhs.num_sulfur != rhs.num_sulfur) return lhs.num_sulfur < rhs.num_sulfur;
      if (lhs.num_comp_sulfur != rhs.num_comp_sulfur) return lhs.num_comp_sulfur < rhs.num_comp_sulfur;
      if (lhs.num_selenium != rhs.num_selenium) return lhs.num_selenium < rhs.num_selenium;
      return lhs.num_comp_selenium < rhs.num_comp_selenium;
    }

    ModelAttributes TensorProductSpline::get_model_attributes() {
      return ModelAttributes(num_sulfur, num_comp_sulfur, num_selenium, num_comp_selenium, precursor_isotope, fragment_isotope);
    }
}
