/** \file ProfileFitter.cxx
    \brief Implementation of ProfileFitter class.
    \author James Peachey, HEASARC/GSSC, based on pulsefitter6.pro by Jay Norris.
*/
#include <algorithm>
#include <stdexcept>

#include <cassert>

#include "burstFit/ProfileFitter.h"

#include "evtbin/Binner.h"
#include "evtbin/Hist1D.h"

#include <iostream>

namespace burstFit {

  const double ProfileFitter::s_fract_threshold = .05;

  ProfileFitter::ProfileFitter(const evtbin::Hist1D * hist): m_peak_index(), m_valley_index(), m_init_fit_par(),
    m_fit_par(), m_hist(hist) {
    findPeaks();
    guessInitialParameters();
  }

  ProfileFitter::~ProfileFitter() {}

  void ProfileFitter::fit() {
  }

  const ProfileFitter::IndexCont_t & ProfileFitter::getPeaks() { return m_peak_index; }

  ProfileFitter::FitParCont_t & ProfileFitter::getInitialParameters() { return m_init_fit_par; }

  const ProfileFitter::FitParCont_t & ProfileFitter::getInitialParameters() const { return m_init_fit_par; }

  const ProfileFitter::FitParCont_t & ProfileFitter::getParameters() const { return m_fit_par; }

  void ProfileFitter::findPeaks() {
    m_peak_index.clear();
    m_valley_index.clear();
    // The following is based on findpvs3.pro:
    if (!m_hist->getBinners().empty() && m_hist->getBinners().front()->getNumBins() >= 3) {
      // Need number of bins.
      Index_t num_bins = m_hist->getBinners().front()->getNumBins();
    
      // Get maximum bin value, which plays a role in determining whether threshold is exceeded.
      double max_bin_val = *std::max_element(m_hist->begin(), m_hist->end());

      double first_bin_val = (*m_hist)[0];
      double last_bin_val = (*m_hist)[num_bins - 1];

      // First and last bin are used to set a background level, which in turn helps set thresholds for determining
      // significance of peaks and valleys.
      double background = .5 * (first_bin_val + last_bin_val);

      // If a maximum bin occurs in the first or last position in the profile, the background is as high
      // or higher than any peak. The background estimates thus would completely wreck this algorithm.
      if (first_bin_val == max_bin_val || last_bin_val == max_bin_val)
        throw std::runtime_error("ProfileFitter::findPeaks(): Background is as large or larger than any peaks");

      // First look for all peaks, which are by definition any/all local maxima.
      // Container to hold the indices of CANDIDATE peak locations.
      IndexCont_t peak_index;
      double total_peak = 0.;
      for (Index_t index = 1; index != num_bins - 1; ++index) {
        if ((*m_hist)[index - 1] < (*m_hist)[index] && (*m_hist)[index] > (*m_hist)[index + 1]) {
          // All local maxima count as peaks.
          peak_index.push_back(index);
  
          // Compute the sum of all peaks for use below.
          total_peak += (*m_hist)[index];
        }
      }

      // Make sure there is at least one peak.
      if (peak_index.empty()) throw std::runtime_error("ProfileFitter::findPeaks(): No peaks were found");

      // Use these peak candidates to find valleys, which in turn will determine the final peak choices.
      // Container to hold the indices of valley locations. Start off with one "virtual" valley which is
      // guaranteed to be before the first peak candidate.
      m_valley_index.reserve(peak_index.size() + 1);

      // Determine the average peak height for use in setting the threshold for valleys.
      double average_peak = total_peak / peak_index.size() - background;

      // Threshold for accepting a virtual valley is a fraction of the average peak height.
      double accept_threshold = s_fract_threshold * average_peak;

      // Find first rise point, that is the bin at which a threshold of the first peak is exceeded.
      Index_t threshold_index = 1;
      for (; threshold_index != peak_index.front() && ((*m_hist)[threshold_index] - first_bin_val) < accept_threshold;
        ++threshold_index) {}

      // First "virtual" valley is the bin before the first rise bin.
      m_valley_index.push_back(threshold_index - 1);

      // Threshold for retaining valleys is a larger fraction of the highest peak height.
      double retain_threshold = .9 * (max_bin_val - background);

      // Find interjacent valleys, which are minima between the peak candidates.
      for (IndexCont_t::iterator itor = peak_index.begin(); itor != peak_index.end() - 1; ++itor) {
        // Find position of the minimum.
        DataCont_t::const_iterator min_pos = std::min_element(m_hist->begin() + *itor, m_hist->begin() + *(itor + 1));

        // Treat minimum as a valley only if it is low enough.
        if (*min_pos - background < retain_threshold) m_valley_index.push_back(min_pos - m_hist->begin());
      }

      // Find final decay point, that is the first bin after the last peak in which the bin value falls below the threshold.
      for (threshold_index = peak_index.back() + 1;
        threshold_index != num_bins && (*m_hist)[threshold_index] - last_bin_val >= accept_threshold; ++threshold_index) {}

      // Final "virtual" valley is the final decay bin.
      m_valley_index.push_back(threshold_index);

      // Final choices for peaks are absolute maxima between significant valleys.
      for (IndexCont_t::iterator itor = m_valley_index.begin(); itor != m_valley_index.end() - 1; ++itor) {
        m_peak_index.push_back(std::max_element(m_hist->begin() + *itor, m_hist->begin() + *(itor + 1)) - m_hist->begin());
      }
    } else {
      throw std::runtime_error("ProfileFitter::findPeaks(): Not enough blocks to find peaks");
    }
  }

  void ProfileFitter::guessInitialParameters() {
    if (!m_hist->getBinners().empty() && m_hist->getBinners().front()->getNumBins() >= 3) {
      // Resize the vector of initial guesses to have the correct number of parameters.
      m_init_fit_par.resize(m_peak_index.size(), FitPar_t(NumParam, 0.));

      // For now, take background to be the average of the first and last bin, just
      // like in findPeaks.
      // TODO: How should the background really be handled?
      const evtbin::Binner * binner = m_hist->getBinners().front();
      size_t num_bins = binner->getNumBins();
      double first_bin_val = (*m_hist)[0];
      double last_bin_val = (*m_hist)[num_bins - 1];
      double background = .5 * (first_bin_val + last_bin_val);
      for (FitParCont_t::size_type par_index = 0; par_index != m_init_fit_par.size(); ++par_index) {
        // Some local aliases to improve readability.
        FitPar_t & fit_par(m_init_fit_par[par_index]);
        Index_t peak_index = m_peak_index[par_index];
        Index_t valley_index = m_valley_index[par_index];

// TODO: This should always be the case, so this assert should not really be necessary.
assert(valley_index < peak_index);

        // Guessed amplitude is just the height above background.
        fit_par[Amplitude] = (*m_hist)[peak_index] - background;
        // Guessed peak center is the point between valley and peak at which the
        // block value rises above the first block's value.
        for (IndexCont_t::size_type index = valley_index; index != peak_index; ++index) {
          if ((*m_hist)[index] > first_bin_val) {
            fit_par[Time0] = binner->getInterval(index).midpoint();
            break;
          }
        }

        // Guessed time constants are the time difference between Time0 and the valley.
        fit_par[Tau1] = binner->getInterval(peak_index).midpoint() - fit_par[Time0];
        fit_par[Tau2] = fit_par[Tau1];

        // Guessed background is the same background used throughout.
        fit_par[Background] = background;
std::clog << "fit_par[Amplitude] is " << fit_par[Amplitude] << std::endl;
std::clog << "fit_par[Time0] is " << fit_par[Time0] << std::endl;
std::clog << "fit_par[Tau1] is " << fit_par[Tau1] << std::endl;
std::clog << "fit_par[Tau2] is " << fit_par[Tau2] << std::endl;
std::clog << "fit_par[Background] is " << fit_par[Background] << std::endl;
      }
    }
  }

}
