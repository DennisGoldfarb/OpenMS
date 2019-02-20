// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#pragma once

// Interfaces
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <OpenMS/KERNEL/FeatureMap.h>

#include <fstream>

namespace OpenMS
{

  /**
   * @brief Class to write out an OpenSwath OSW SQLite output (PyProphet input)
   *
   * The class can take a FeatureMap and create a set of string from it
   * suitable for output to OSW using the prepareLine function.
   *
   */
  class OPENMS_DLLAPI OpenSwathOSWWriter
  {
    String output_filename_;
    String input_filename_;
    OpenMS::UInt64 run_id_;
    bool doWrite_;
    bool use_ms1_traces_;
    bool sonar_;
    bool enable_uis_scoring_;

  public:

    OpenSwathOSWWriter(const String& output_filename,
                       const String& input_filename = "inputfile",
                       bool ms1_scores = false,
                       bool sonar = false,
                       bool uis_scores = false) :
      output_filename_(output_filename),
      input_filename_(input_filename),
      run_id_(OpenMS::UniqueIdGenerator::getUniqueId()),
      doWrite_(!output_filename.empty()),
      use_ms1_traces_(ms1_scores),
      sonar_(sonar),
      enable_uis_scoring_(uis_scores)
      {}

    bool isActive() const;

    /**
     * @brief Initializes file by generating SQLite tables
     *
     */
    void writeHeader();

    /**
     * @brief Prepare scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param feature The feature being evaluated
     * @param score_name The name of the queried score
     *
     * @returns A string with the queried score
     *
     */
    String getScore(const Feature& feature, std::string score_name) const;

    /**
     * @brief Prepare concatenated scores for SQLite insertion
     *
     * Some scores might not be defined, those are reported as NULL
     *
     * @param feature The feature being evaluated
     * @param score_name The name of the queried score
     *
     * @returns A vector of strings with the queried scores
     *
     */
    std::vector<String> getSeparateScore(const Feature& feature, std::string score_name) const;

    /**
     * @brief Prepare a single line (feature) for output
     *
     * The result can be flushed to disk using writeLines (either line by line
     * or after collecting several lines).
     *
     * @param pep The compound (peptide/metabolite) used for extraction
     * @param transition The transition used for extraction 
     * @param output The feature map containing all features (each feature will generate one entry in the output)
     * @param id The transition group identifier (peptide/metabolite id)
     *
     * @returns A string to be written using writeLines
     *
     */
    String prepareLine(const OpenSwath::LightCompound& /* pep */,
        const OpenSwath::LightTransition* /* transition */,
        FeatureMap& output, String id) const;

    /**
     * @brief Write data to disk
     *
     * Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
     * 
     * @param to_osw_output Statements generated by prepareLine
     *
     * @note Try to call this function as little as possible (it opens a new
     * database connection each time)
     *
     * @note Only call inside an OpenMP critical section
     *
     */
    void writeLines(const std::vector<String>& to_osw_output);

  };

}

