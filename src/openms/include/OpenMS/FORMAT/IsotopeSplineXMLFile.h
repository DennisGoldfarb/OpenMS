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

#ifndef OPENMS_FORMAT_ISOTOPESPLINEXMLFILE_H
#define OPENMS_FORMAT_ISOTOPESPLINEXMLFILE_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>


#include <vector>

namespace OpenMS
{
    /**
      @brief Used to load (storing not supported, yet) ProtXML files
      This class is used to load (storing not supported, yet) documents that implement
      the schema of ProtXML files.
          A documented schema for this format comes with the TPP.
          OpenMS can only read parts of the protein_summary subtree to extract protein-peptide associations. All other parts are silently ignored.
          For protein groups, only the "group leader" (which is annotated with a probability and coverage)
          receives these attributes. All indistinguishable proteins of the same group only have
          an accession and score of -1.
      @note This format will eventually be replaced by the HUPO-PSI (mzIdentML and mzQuantML) AnalysisXML formats!
      @todo Document which metavalues of Protein/PeptideHit are filled when reading ProtXML (Chris)
      @ingroup FileIO
    */
    class OPENMS_DLLAPI IsotopeSplineXMLFile :
            protected Internal::XMLHandler,
            public Internal::XMLFile {
    public:

        /// Constructor
        IsotopeSplineXMLFile();

        /**
            @brief Loads the identifications of an ProtXML file without identifier
            The information is read in and the information is stored in the
            corresponding variables
            @exception Exception::FileNotFound is thrown if the file could not be opened
            @exception Exception::ParseError is thrown if an error occurs during parsing
        */
        void load(const String &filename, std::vector<CubicSpline2d>* models);

    protected:

        /// reset members after reading/writing
        void resetMembers_();

        /// Docu in base class
        void endElement(const XMLCh *const /*uri*/, const XMLCh *const /*local_name*/, const XMLCh *const qname);

        /// Docu in base class
        void startElement(const XMLCh *const /*uri*/, const XMLCh *const /*local_name*/, const XMLCh *const qname,
                          const xercesc::Attributes &attributes);

        /// Docu in base class
        void characters(const XMLCh * const chars, const XMLSize_t length);


        Base64 decoder_;

        /// @name members for loading data
        //@{
        /// model
        std::vector<CubicSpline2d> *models_;
        UInt num_models_;
        Int num_sulfur_;
        UInt isotope_;
        Base64::ByteOrder byteOrder_;
        UInt precision_;
        UInt length_;
        UInt spline_order_;
        String base64Data_;
        std::vector<double> a_;
        std::vector<double> b_;
        std::vector<double> c_;
        std::vector<double> d_;
        std::vector<double> x_;
        //@}

    private:
        void processBinaryArrayAttributes_(const xercesc::Attributes &attributes);

    };
} // namespace OpenMS

#endif //OPENMS_FORMAT_ISOTOPESPLINEXMLFILE_H
