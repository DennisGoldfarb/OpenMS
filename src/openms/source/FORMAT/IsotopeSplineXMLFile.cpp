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


#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/IsotopeSplineXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

namespace OpenMS
{

    IsotopeSplineXMLFile::IsotopeSplineXMLFile() :
            XMLHandler("", "1.2"),
            XMLFile("/SCHEMAS/IsotopeSplineXML_1_0.xsd", "1.0")
    {
    }

    void IsotopeSplineXMLFile::load(const String& filename, std::map<ModelAttributes, TensorProductSpline*>* models)
    {
      //Filename for error messages in XMLHandler
      file_ = filename;
      this->models = models;

      resetMembers_();

      parse_(filename, this);
    }

    /// reset members
    void IsotopeSplineXMLFile::resetMembers_()
    {

    }

    void IsotopeSplineXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
    {
      String tag = sm_.convert(qname);
      open_tags_.push_back(tag);

      if (tag == "model")
      {
        model = new TensorProductSpline();
        model->num_sulfur = 0;
        model->num_comp_sulfur = 0;
        model->num_selenium = 0;
        model->num_comp_selenium = 0;
        optionalAttributeAsInt_(model->num_sulfur, attributes, "S");
        optionalAttributeAsInt_(model->num_comp_sulfur, attributes, "CS");
        optionalAttributeAsInt_(model->num_selenium, attributes, "Se");
        optionalAttributeAsInt_(model->num_comp_selenium, attributes, "CSe");

        model->precursor_isotope = attributeAsInt_(attributes, "PrecursorIsotope");
        model->fragment_isotope = attributeAsInt_(attributes, "FragmentIsotope");

        model->order_precursor_mass = attributeAsInt_(attributes, "PrecursorOrder");
        model->order_fragment_mass = attributeAsInt_(attributes, "FragmentOrder");
      }
      else if (tag == "fragmentMassBreaks" || tag == "precursorMassBreaks" || tag == "coefficients")
      {
        processBinaryArrayAttributes_(attributes);
      }
    }

    void IsotopeSplineXMLFile::characters(const XMLCh * const chars, const XMLSize_t length)
    {
      if (open_tags_.back() == "fragmentMassBreaks")
      {
        // Since we convert a Base64 string here, it can only contain plain ASCII
        sm_.appendASCII(chars, length, base64Data);
        decoder_.decode(base64Data, byteOrder, model->breaks_fragment_mass);
      }
      else if (open_tags_.back() == "precursorMassBreaks")
      {
        // Since we convert a Base64 string here, it can only contain plain ASCII
        sm_.appendASCII(chars, length, base64Data);
        decoder_.decode(base64Data, byteOrder, model->breaks_precursor_mass);
      }
      else if (open_tags_.back() == "coefficients")
      {
        std::vector<float> data;
        sm_.appendASCII(chars, length, base64Data);
        decoder_.decode(base64Data, byteOrder, data);

        int num_coefficients = model->order_fragment_mass * model->order_precursor_mass;

        // allocate memory
        UInt num_breaks_precursor_mass = model->breaks_precursor_mass.size();
        UInt num_breaks_fragment_mass = model->breaks_fragment_mass.size();

        model->coefficients = new float**[num_breaks_precursor_mass];
        int max_j = 0;
        for (int i = 0; i < num_breaks_precursor_mass; ++i) {
          for (; max_j < num_breaks_fragment_mass &&  model->breaks_fragment_mass[max_j] <  model->breaks_precursor_mass[i]; ++max_j) {}
          model->coefficients[i] = new float*[max_j+1];
          for (int j = 0; j <= max_j; ++j) {
            model->coefficients[i][j] = new float[num_coefficients];
          }
        }


        int p_index=0, f_index=-1, c_index=0, last_f_index=0;
        float precursor_mass =  model->breaks_precursor_mass[p_index], fragment_mass;

        for (int i = 0; i < data.size(); i++) {
          c_index = i%num_coefficients;

          if (c_index==0) {
            f_index++;
            fragment_mass =  model->breaks_fragment_mass[f_index];

            if (fragment_mass > precursor_mass || f_index == num_breaks_fragment_mass) {
              f_index=0;
              p_index++;
              precursor_mass =  model->breaks_precursor_mass[p_index];
            }
          }

          model->coefficients[p_index][f_index][c_index] = data[i];
        }
      }
    }

    void IsotopeSplineXMLFile::processBinaryArrayAttributes_(const xercesc::Attributes& attributes)
    {
      byteOrder = attributeAsString_(attributes, "endian") == String("little") ? Base64::BYTEORDER_LITTLEENDIAN : Base64::BYTEORDER_BIGENDIAN;
      precision = attributeAsInt_(attributes, "precision");
      length = attributeAsInt_(attributes, "length");
    }

    void IsotopeSplineXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
    {
      String tag = sm_.convert(qname);
      open_tags_.pop_back();

      if (tag == "model")
      {
        models->insert(std::make_pair(model->get_model_attributes(), model));
      }
      else if (tag == "fragmentMassBreaks" || tag == "precursorMassBreaks" || tag == "coefficients")
      {
        base64Data.clear();
      }
    }


} // namespace OpenMS

