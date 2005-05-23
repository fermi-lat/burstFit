/** \file ProfileFitter.h
    \brief Declaration of ProfileFitter class.
*/
#ifndef burstFit_ProfileFitter_h
#define burstFit_ProfileFitter_h

#include <vector>

namespace evtbin {
  class Hist1D;
}

namespace burstFit {

  /** \class ProfileFitter
      \brief Class for fitting a binned data profile to a function.
  */
  class ProfileFitter {
    public:
      typedef std::vector<double> DataCont_t;
      typedef DataCont_t::size_type Index_t;
      typedef std::vector<Index_t> IndexCont_t;
      typedef std::vector<double> FitPar_t;
      typedef std::vector<FitPar_t> FitParCont_t;

      ProfileFitter(const evtbin::Hist1D * hist);

      virtual ~ProfileFitter();

      virtual void fit();

      virtual const IndexCont_t & getPeaks();

      virtual FitParCont_t & getInitialParameters();

      virtual const FitParCont_t & getInitialParameters() const;

      virtual const FitParCont_t & getParameters() const;

    protected:
      enum ParamId {
        Amplitude,
        Time0,
        Tau1,
        Tau2,
        Background,
        NumParam
      };

      static const double s_fract_threshold;
      virtual void findPeaks();

      virtual void guessInitialParameters();

    private:
      IndexCont_t m_peak_index;
      IndexCont_t m_valley_index;
      FitParCont_t m_init_fit_par;
      FitParCont_t m_fit_par;
      const evtbin::Hist1D * m_hist;
  };

}

#endif
