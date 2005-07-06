/** \file BurstModel.h
    \brief Implementation of BurstModel class.
*/
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "burstFit/BurstModel.h"

#include "evtbin/Binner.h"
#include "evtbin/Hist1D.h"

#include "optimizers/Parameter.h"
#include "optimizers/dArg.h"

#include "st_stream/Stream.h"

namespace {
  enum { Amplitude, Time0, Tau1, Tau2, NumParPerPeak };
}

namespace burstFit {

  const double BurstModel::s_fract_threshold = .05;

  BurstModel::BurstModel(const FitPar_t & parameter): m_peak_index(), m_valley_index() {
    if (1 != (parameter.size() % 4))
      throw std::logic_error("BurstModel::BurstModel(parameter): There must be 4 parameters per peak, plus one background term.");

    setMaxNumParams(parameter.size());
    setParams(const_cast<FitPar_t &>(parameter));

    m_funcType = Addend;
    m_argType = "dArg";
    m_genericName = "BurstModel";
  }

  BurstModel::BurstModel(const evtbin::Hist1D * hist): m_peak_index(), m_valley_index() {
    findPeaks(hist);
    guessInitialParameters(hist, m_parameter);
    setMaxNumParams(m_parameter.size());

    m_funcType = Addend;
    m_argType = "dArg";
    m_genericName = "BurstModel";
  }

  double BurstModel::value(optimizers::Arg & x) const {
    using namespace optimizers;
    double abscissa = dynamic_cast<dArg &>(x).getValue();

    double value = 0.;
    for (size_t index = 0; index != m_parameter.size() / 4; ++index) {
      // The form of this computation is taken from dofit_prodexpback.pro.
      double t = abscissa - m_parameter[4 * index + Time0].getTrueValue();
      if (t <= 0.) continue;
      double amplitude = m_parameter[4 * index + Amplitude].getTrueValue();
      double coeff1 = m_parameter[4 * index + Tau1].getTrueValue();
      double coeff2 = m_parameter[4 * index + Tau2].getTrueValue();

      double mu = sqrt(coeff1 / coeff2);
      // In findPVs, lambda was used to help compute parameters. Do not use lambda here, because 2 * mu can
      // be large enough to blow up the exponential even if 2 * mu - B is not.
      // double lambda = exp(2. * mu);
      // value += amplitude * lambda * exp(-B);
      double B = coeff1 / t + t / coeff2;
      value += amplitude * exp(2. * mu - B);
    }
    // Add the flat background, which is the last parameter.
    value += m_parameter.back().getTrueValue();

    return value;
  }

  double BurstModel::derivByParam(optimizers::Arg & x, const std::string & par_name) const {
    using namespace optimizers;
    if (std::string::npos != par_name.find("Bckgnd")) return 1.;
    double abscissa = dynamic_cast<dArg &>(x).getValue();
    // Do this just to make sure the parameter name given is valid. getParam will throw if
    // it is not.
    getParam(par_name);

    // Pull out the number from the parameter name.
    std::stringstream ss;
    ss << par_name.substr(4, par_name.find("::"));

    size_t index;
    ss >> index;

    double t = abscissa - m_parameter[4 * index + Time0].getTrueValue();
    if (t <= 0.) return 0.;
    double amplitude = m_parameter[4 * index + Amplitude].getTrueValue();
    double coeff1 = m_parameter[4 * index + Tau1].getTrueValue();
    double coeff2 = m_parameter[4 * index + Tau2].getTrueValue();

    double mu = sqrt(coeff1 / coeff2);
    // In findPVs, lambda was used to help compute parameters. Do not use lambda here, because 2 * mu can
    // be large enough to blow up the exponential even if 2 * mu - B is not.
    // double lambda = exp(2. * mu);
    // double factor = amplitude * lambda * exp(-B);
    double B = coeff1 / t + t / coeff2;

    double factor = amplitude * exp(2. * mu - B);

    double partial = 0.;
    if (std::string::npos != par_name.find("Amp")) partial = factor / amplitude;
    else if (std::string::npos != par_name.find("Time0")) partial = factor * (1. / coeff2 - coeff1 / (t * t));
    else if (std::string::npos != par_name.find("Tau1")) partial = factor * (1. / sqrt(coeff1 * coeff2) - 1. / t);
    else if (std::string::npos != par_name.find("Tau2")) partial = factor * (t / (coeff2 * coeff2) - mu / coeff2);
    return partial;
  }

  optimizers::Function * BurstModel::clone() const { return new BurstModel(*this); }

  st_stream::OStream & BurstModel::write(st_stream::OStream & os) const {
    for (unsigned int index = 0; index != m_parameter.size() / 4; ++index) {
      os << "Peak " << index + 1 << ":" << std::endl;
      //os << std::right();
      os.width(17); os << "Amplitude = " << m_parameter[4 * index + Amplitude].getTrueValue() << std::endl;
      os.width(17); os << "Time0 = " << m_parameter[4 * index + Time0].getTrueValue() << std::endl;
      os.width(17); os << "Tau1 = " << m_parameter[4 * index + Tau1].getTrueValue() << std::endl;
      os.width(17); os << "Tau2 = " << m_parameter[4 * index + Tau2].getTrueValue() << std::endl;
    }
    os.width(17); os << "Background = " << m_parameter.back().getTrueValue();
    return os;
  }

  void BurstModel::findPeaks(const evtbin::Hist1D * hist) {
    m_peak_index.clear();
    m_valley_index.clear();
    // The following is based on findpvs3.pro:
    if (!hist->getBinners().empty() && hist->getBinners().front()->getNumBins() >= 3) {
      // Need number of bins.
      Index_t num_bins = hist->getBinners().front()->getNumBins();
    
      // Get maximum bin value, which plays a role in determining whether threshold is exceeded.
      double max_bin_val = *std::max_element(hist->begin(), hist->end());

      double first_bin_val = (*hist)[0];
      double last_bin_val = (*hist)[num_bins - 1];

      // First and last bin are used to set a background level, which in turn helps set thresholds for determining
      // significance of peaks and valleys.
      double background = .5 * (first_bin_val + last_bin_val);

      // If a maximum bin occurs in the first or last position in the profile, the background is as high
      // or higher than any peak. The background estimates thus would completely wreck this algorithm.
      if (first_bin_val == max_bin_val || last_bin_val == max_bin_val)
        throw std::runtime_error("BurstModel::findPeaks(): Background is as large or larger than any peaks");

      // First look for all peaks, which are by definition any/all local maxima.
      // Container to hold the indices of CANDIDATE peak locations.
      IndexCont_t peak_index;
      double total_peak = 0.;
      for (Index_t index = 1; index != num_bins - 1; ++index) {
        if ((*hist)[index - 1] < (*hist)[index] && (*hist)[index] > (*hist)[index + 1]) {
          // All local maxima count as peaks.
          peak_index.push_back(index);
  
          // Compute the sum of all peaks for use below.
          total_peak += (*hist)[index];
        }
      }

      // Make sure there is at least one peak.
      if (peak_index.empty()) throw std::runtime_error("BurstModel::findPeaks(): No peaks were found");

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
      for (; threshold_index != peak_index.front() && ((*hist)[threshold_index] - first_bin_val) < accept_threshold;
        ++threshold_index) {}

      // First "virtual" valley is the bin before the first rise bin.
      m_valley_index.push_back(threshold_index - 1);

      // Threshold for retaining valleys is a larger fraction of the highest peak height.
      double retain_threshold = .9 * (max_bin_val - background);

      // Find interjacent valleys, which are minima between the peak candidates.
      for (IndexCont_t::iterator itor = peak_index.begin(); itor != peak_index.end() - 1; ++itor) {
        // Find position of the minimum.
        DataCont_t::const_iterator min_pos = std::min_element(hist->begin() + *itor, hist->begin() + *(itor + 1));

        // Treat minimum as a valley only if it is low enough.
        if (*min_pos - background < retain_threshold) m_valley_index.push_back(min_pos - hist->begin());
      }

      // Find final decay point, that is the first bin after the last peak in which the bin value falls below the threshold.
      for (threshold_index = peak_index.back() + 1;
        threshold_index != num_bins && (*hist)[threshold_index] - last_bin_val >= accept_threshold; ++threshold_index) {}

      // Final "virtual" valley is the final decay bin.
      m_valley_index.push_back(threshold_index);

      // Final choices for peaks are absolute maxima between significant valleys.
      for (IndexCont_t::iterator itor = m_valley_index.begin(); itor != m_valley_index.end() - 1; ++itor) {
        m_peak_index.push_back(std::max_element(hist->begin() + *itor, hist->begin() + *(itor + 1)) - hist->begin());
      }
    } else {
      throw std::runtime_error("BurstModel::findPeaks(): Not enough blocks to find peaks");
    }
  }

  void BurstModel::guessInitialParameters(const evtbin::Hist1D * hist, FitPar_t & parameter) const {
    using namespace optimizers;
    parameter.clear();
    if (!hist->getBinners().empty() && hist->getBinners().front()->getNumBins() >= 3) {
      // For now, take background to be the average of the first and last bin, just
      // like in findPeaks.
      const evtbin::Binner * binner = hist->getBinners().front();
      Index_t num_bins = binner->getNumBins();
      double first_bin_val = (*hist)[0];
      double last_bin_val = (*hist)[num_bins - 1];
      double background = .5 * (first_bin_val + last_bin_val);

      parameter.resize(m_peak_index.size() * NumParPerPeak + 1);
      for (IndexCont_t::size_type term_index = 0; term_index != m_peak_index.size(); ++term_index) {
        // Some local aliases to improve readability.
        FitPar_t::size_type par_index = term_index * 4;
        Index_t peak_index = m_peak_index[term_index];
        Index_t valley_index = m_valley_index[term_index];

        // String used to name parameters for this peak.
        std::ostringstream os;
        os << "_" << term_index;

        // Guessed amplitude is just the height above background.
        double amplitude = (*hist)[peak_index] - background;
        parameter[par_index + Amplitude] = Parameter("Amp" + os.str(), amplitude, true);
        parameter[par_index + Amplitude].setBounds(0., 2. * amplitude);

        // Guessed peak center is the point between valley and peak at which the
        // block value rises above the first block's value.
        IndexCont_t::size_type center_index = valley_index;
        for (; center_index < peak_index - 1 && (*hist)[center_index] <= first_bin_val; ++center_index) {}
        double origin = binner->getInterval(center_index).midpoint();
        parameter[par_index + Time0] = Parameter("Time0" + os.str(), origin, true);
        parameter[par_index + Time0].setBounds(binner->getInterval(valley_index).begin(), binner->getInterval(peak_index).end());

        // Guessed time constants are the time differences between peak and origin.
        double coeff = binner->getInterval(peak_index).midpoint() - origin;
        parameter[par_index + Tau1] = Parameter("Tau1" + os.str(), coeff, true);
        parameter[par_index + Tau1].setBounds(0., 3. * coeff);
        parameter[par_index + Tau2] = Parameter("Tau2" + os.str(), coeff, true);
        parameter[par_index + Tau2].setBounds(0., 3. * coeff);

      }
      // Guessed background is the same background used throughout.
      parameter.back() = Parameter("Bckgnd", background, true);
      parameter.back().setBounds(0., 3. * background);

    }
  }

  st_stream::OStream & operator <<(st_stream::OStream & os, const BurstModel & model) {
    return model.write(os);
  }
}
