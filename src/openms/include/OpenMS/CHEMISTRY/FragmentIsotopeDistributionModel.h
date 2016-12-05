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

#ifndef OPENMS_CHEMISTRY_FRAGMENTISOTOPEDISTRIBUTIONMODEL_H
#define OPENMS_CHEMISTRY_FRAGMENTISOTOPEDISTRIBUTIONMODEL_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>

#include <OpenMS/config.h>

namespace OpenMS
{
    /** @ingroup Chemistry

      @brief holds and evaluates tensor product spline models for peptide fragment isotope distributions

      The models stored in this DB are defined in a
      XML file under data/CHEMISTRY/peptideFragmentIsotopeSplines.xml
    */
    class OPENMS_DLLAPI FragmentIsotopeDistributionModel {

    public:

        /// this member function serves as a replacement of the constructor
        inline static FragmentIsotopeDistributionModel* getInstance()
        {
            static FragmentIsotopeDistributionModel* model_ = 0;
            if (model_ == 0)
            {
                model_ = new FragmentIsotopeDistributionModel;
            }
            return model_;
        }

        /** @name Constructors and Destructors
        */
        //@{
        /// destructor
        virtual ~FragmentIsotopeDistributionModel();
        //@}

        double getProbabilities(IsotopeDistribution::ContainerType result, double average_weight_precursor, double average_weight_fragment, const std::vector<UInt>& precursor_isotopes);

        bool canHandleRequest(double average_weight_precursor, double average_weight_fragment, const std::vector<UInt>& precursor_isotopes);
    private:

    protected:

        /** @name Private Constructors
        */
        //@{
        /// default constructor
        FragmentIsotopeDistributionModel();

        ///copy constructor
        FragmentIsotopeDistributionModel(const FragmentIsotopeDistributionModel& fragment_isotope_model);
        //@}

        /** @name Assignment
        */
        //@{
        /// assignment operator
        FragmentIsotopeDistributionModel& operator=(const FragmentIsotopeDistributionModel& fragment_isotope_model);
        //@}

        /**
        @brief reads spline models from the given file

        @throw Exception::ParseError if the file cannot be parsed
        */
        void readSplineModelsFromFile_(const String& filename);

        /// deletes all sub-instances of the stored data
        void clear_();

        /// deletes all models
        void clearModels_();
    };

}

#endif //OPENMS_CHEMISTRY_FRAGMENTISOTOPEDISTRIBUTIONMODEL_H
