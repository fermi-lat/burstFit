#include "burstFit/NegativeStat.h"

namespace burstFit {

  NegativeStat::NegativeStat(optimizers::Statistic * stat): m_stat(stat) {
    // Mirror the supplied statistic's parameters in this object's parameters, so that methods will work correctly
    // when called for this object.
    m_stat->getParams(m_parameter);
  }

  double NegativeStat::value() const {
    return -m_stat->value();
  }

  void NegativeStat::getFreeDerivs(std::vector<double> & derivs) const {
    m_stat->getFreeDerivs(derivs);
    for (std::vector<double>::iterator itor = derivs.begin(); itor != derivs.end(); ++itor) *itor = -*itor;
  }

  std::vector<double>::const_iterator NegativeStat::setFreeParamValues_(std::vector<double>::const_iterator it) {
    m_stat->setFreeParamValues_(it);
    // Also store the parameters in this object's parameter set, so that other methods from Function such
    // as getFreeParamValues() work correctly .
    return Function::setFreeParamValues_(it);
  }

  optimizers::Function * NegativeStat::clone() const { return new NegativeStat(m_stat); }

  double NegativeStat::value(optimizers::Arg & x) const {
    return -dynamic_cast<const optimizers::Function *>(m_stat)->value(x);
  }

  double NegativeStat::derivByParam(optimizers::Arg & x, const std::string & parameter_name) const {
    return -m_stat->derivByParam(x, parameter_name);
  }

  void NegativeStat::getFreeDerivs(optimizers::Arg & x, std::vector<double> & derivs) const {
    dynamic_cast<const optimizers::Function *>(m_stat)->getFreeDerivs(x, derivs);
    for (std::vector<double>::iterator itor = derivs.begin(); itor != derivs.end(); ++itor) *itor = -*itor;
  }
}
