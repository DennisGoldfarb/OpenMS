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


#include <OpenMS/FORMAT/IsotopeSplineXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{

    IsotopeSplineXMLFile::IsotopeSplineXMLFile() :
            XMLHandler("", "1.2"),
            XMLFile("/SCHEMAS/IsotopeSplineXML_1_0.xsd", "1.0")
    {
    }

    void IsotopeSplineXMLFile::load(const String& filename)
    {
      //Filename for error messages in XMLHandler
      file_ = filename;

      resetMembers_();

      parse_(filename, this);
    }

    /// reset members
    void IsotopeSplineXMLFile::resetMembers_()
    {
      a_.clear();
      b_.clear();
      c_.clear();
      d_.clear();
      x_.clear();
      base64Data_.clear();
    }

    void IsotopeSplineXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      String tag = sm_.convert(qname);
      open_tags_.push_back(tag);
      if (tag == "models")
      {
        max_isotope_ = attributeAsInt_(attributes, "maxIsotope");
        max_sulfur_ = attributeAsInt_(attributes, "maxSulfur");
        models_.resize(max_isotope_ + 1);
        sulfur_specific_models_.resize((max_isotope_ + 1) * (max_sulfur_ + 1));
      }
      else if (tag == "model")
      {
        is_sulfur_model_ = optionalAttributeAsUInt_(num_sulfur_, attributes, "S");
        isotope_ = attributeAsInt_(attributes, "isotope");
        spline_order_ = attributeAsInt_(attributes, "order");
        if (spline_order_ != 4)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Only splines of order 4 are currently supported.");
        }
      }
      else if (tag == "knots" || tag == "coefficients")
      {
        processBinaryArrayAttributes_(attributes);
      }
    }

    void IsotopeSplineXMLFile::characters(const XMLCh * const chars, const XMLSize_t length)
    {
      if (open_tags_.back() == "knots")
      {
        // Since we convert a Base64 string here, it can only contain plain ASCII
        sm_.appendASCII(chars, length, base64Data_);
        decoder_.decode(base64Data_, byteOrder_, x_);
      }
      else if (open_tags_.back() == "coefficients")
      {
        std::vector<double> data;
        sm_.appendASCII(chars, length, base64Data_);
        decoder_.decode(base64Data_, byteOrder_, data);

        for (UInt i = 0; i < data.size(); i+=spline_order_)
        {
          a_.push_back(data[i]);
          b_.push_back(data[i+1]);
          c_.push_back(data[i+2]);
          d_.push_back(data[i+3]);
        }
      }
    }

    void IsotopeSplineXMLFile::processBinaryArrayAttributes_(const xercesc::Attributes& attributes)
    {
      byteOrder_ = attributeAsString_(attributes, "endian") == String("little") ? Base64::BYTEORDER_LITTLEENDIAN : Base64::BYTEORDER_BIGENDIAN;
      precision_ = attributeAsInt_(attributes, "precision");
      length_ = attributeAsInt_(attributes, "length");
    }

    void IsotopeSplineXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      String tag = sm_.convert(qname);
      open_tags_.pop_back();

      if (tag == "model")
      {
        CubicSpline2d spline = CubicSpline2d(a_, b_, c_, d_, x_);
        if (is_sulfur_model_)
        {
          sulfur_specific_models_[isotope_ + (num_sulfur_ * (max_isotope_ + 1))] = spline;
        } else
        {
          models_[isotope_] = spline;
        }

        resetMembers_();
      }
      else if (tag == "knots" || tag == "coefficients")
      {
        base64Data_.clear();
      }
    }

    UInt IsotopeSplineXMLFile::getMaxIsotope() const
    {
      return max_isotope_;
    }

    UInt IsotopeSplineXMLFile::getMaxSulfur() const
    {
      return max_sulfur_;
    }

    std::vector<CubicSpline2d> IsotopeSplineXMLFile::getModels() const
    {
      return models_;
    }

    std::vector<CubicSpline2d> IsotopeSplineXMLFile::getSulfurSpecificModels() const
    {
      return sulfur_specific_models_;
    }

} // namespace OpenMS
