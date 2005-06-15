/** \file BurstModel.h
    \brief Declaration of BurstModel class.
*/
#ifndef burstFit_BurstModel_h
#define burstFit_BurstModel_h

#include <vector>

#include "optimizers/Function.h"

namespace evtbin {
  class Hist1D;
}

namespace optimizers {
  class Arg;
  class Parameter;
}

namespace burstFit {

  /** \class BurstModel
      \brief Function which is a sum of any number of terms of the form:
      A * exp( - ( (x - origin) / coeff1 + coeff2 / (x - origin) ) )
      plus a single constant background term.
  */
  class BurstModel : public optimizers::Function {
    public:
      typedef std::vector<double> DataCont_t;
      typedef DataCont_t::size_type Index_t;
      typedef std::vector<Index_t> IndexCont_t;
      typedef std::vector<double> FitPar_t;

      BurstModel(const FitPar_t & parameter);

      BurstModel(const evtbin::Hist1D * hist);

      virtual double value(optimizers::Arg & x) const;

      virtual double derivByParam(optimizers::Arg & x, const std::string & par_name) const;

      virtual optimizers::Function * clone() const;

    protected:
      static const double s_fract_threshold;

      virtual void findPeaks(const evtbin::Hist1D * hist);

      virtual void guessInitialParameters(const evtbin::Hist1D * hist, FitPar_t & parameter) const;

      virtual void setParameters(const FitPar_t & parameter);

    private:
      IndexCont_t m_peak_index;
      IndexCont_t m_valley_index;
  };

}
#endif
