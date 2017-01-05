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


#ifndef OPENMS_HOST_ISOTOPESPLINEDB_H
#define OPENMS_HOST_ISOTOPESPLINEDB_H

#include <map>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  struct OPENMS_DLLAPI ModelAttributes {
    Int num_sulfur;
    Int isotope;

    ModelAttributes() {}

    ModelAttributes(Int S, Int i) {
      num_sulfur = S;
      isotope = i;
    }
  };

  bool operator<(const ModelAttributes & lhs, const ModelAttributes & rhs);

  /** @ingroup Chemistry
   *  @brief holds and evaluates 2D cubic spline models for peptide isotope distributions
   *  The models stored in this DB are defined in an
   *  XML file under data/CHEMISTRY/IsotopeSplines.xml
   */
  class OPENMS_DLLAPI IsotopeSplineDB
  {

    public:

      typedef std::map<ModelAttributes, CubicSpline2d*>::iterator Iterator;

      /// this member function serves as a replacement of the constructor
      inline static IsotopeSplineDB* getInstance()
      {
        static IsotopeSplineDB* db_ = 0;
        if (db_ == 0)
        {
          db_ = new IsotopeSplineDB;
        }
        return db_;
      }

      /** @name Constructors and Destructors
      */
      //@{
      /// destructor
      virtual ~IsotopeSplineDB();
      //@}

      void approximateIsotopeDistribution(IsotopeDistribution::ContainerType& result, double average_weight, UInt max_isotope);

      bool inModelBounds(double average_weight, UInt max_isotope);
    private:

    protected:

      /** @name Private Constructors
      */
      //@{
      /// default constructor
      IsotopeSplineDB();

      ///copy constructor
      IsotopeSplineDB(const IsotopeSplineDB& isotope_spline_db);
      //@}

      /** @name Assignment
      */
      //@{
      /// assignment operator
      IsotopeSplineDB& operator=(const IsotopeSplineDB& isotope_spline_db);
      //@}

      /**
      @brief reads spline models from the given file
      @throw Exception::ParseError if the file cannot be parsed
      */
      void readSplinesFromFile_(const String& filename);

      /// deletes all sub-instances of the stored data
      void clear_();

      /// deletes all models
      void clearModels_();

      double min_mass_;
      double max_mass_;

      UInt max_isotope_;

      std::map<ModelAttributes, CubicSpline2d*> models;
  };

}


#endif //OPENMS_HOST_ISOTOPESPLINEDB_H
