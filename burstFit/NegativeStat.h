#ifndef burstFit_NegativeStat_h
#define burstFit_NegativeStat_h

#include <vector>

#include "optimizers/Statistic.h"

namespace optimizers {
  class Arg;
  class Function;
}

namespace burstFit {

  class NegativeStat : public optimizers::Statistic {
    public:
      NegativeStat(optimizers::Statistic * stat);

      // From Statistic:
      virtual double value() const;

      // From Statistic:
      virtual void getFreeDerivs(std::vector<double> & derivs) const;

      // From Function:
      virtual std::vector<double>::const_iterator setFreeParamValues_(std::vector<double>::const_iterator it);

      virtual optimizers::Function * clone() const;

    protected:
      // From Function:
      virtual double value(optimizers::Arg & x) const;

      // From Function:
      virtual double derivByParam(optimizers::Arg & x, const std::string & parameter_name) const;

      // From Function:
      virtual void getFreeDerivs(optimizers::Arg & x, std::vector<double> & derivs) const;

    private:
      optimizers::Statistic * m_stat;
  };

}

#endif
